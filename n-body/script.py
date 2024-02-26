import numpy as np
import pandas as pd
import time

import copy # Usado para criar cópias não vinculadas do vetor de estado
from math import ceil # Útil para arredondar floats para inteiros
from itertools import combinations # Produz os conjuntos de combinações para pares de números

class nBody:
    def __init__(self, Si=None, masses=None, eps=0, G=6.67408e-11, folder='', interSnaps=0):
        # Vetor de estado e massas
        self.Si = Si # Define o vetor de estado inicial
        self.masses = masses # Define o array de massas do sistema

        # Define a constante gravitacional
        self.G = G

        # Historico dos valores da simulação
        self.history = None
        self.energies = None

        # Historico dos valores de posição no referêncial do CM
        self.positions = None

        # Variáveis úteis
        self.N, self.D = Si.shape # Número de corpos e dimensionalidade do vetor de estado
        self.D = self.D // 2 # Dimensão do problema, 2D ou 3D

        # Comprimento de amortecimento
        self.eps = eps

        # Diretório final:
        self.dir = folder

        # Número de arquivos de saída:
        self.nSnp = interSnaps

    def get_energy(self, S):
        '''
        Entrada:
        - S: vetor de estado atual com posições e velocidades
        Saída:
        - Energia cinética e potencial do estado
        '''

        R = S[:, :self.D] # Submatriz com todas as posições
        dotR = S[:, self.D:] # Submatriz com todas as velocidades
        masses = copy.deepcopy(self.masses).reshape(self.N,1)

        # Energia Cinética:
        KE = 0.5 * np.sum(np.sum(masses * dotR**2 ))

        # Energia Potencial
        PE = 0
        for body_i, body_j in self.pairs:
            r = np.linalg.norm(R[body_j] - R[body_i]) # Distancia entre os corpos
            PE -= self.masses[body_i] * self.masses[body_j] / r
        PE *= self.G

        return KE, PE

    def get_state_deriv(self, S):
        '''
        Entrada:
        - S: vetor de estado atual com posições e velocidades
        Saída:
        - dotS: vetor derivada de S com velocidades e acelerações
        '''

        R = S[:, :self.D] # Submatriz com todas as posições
        dotR = S[:, self.D:] # Submatriz com todas as velocidades

        # Constrói a estrutura do Output
        dotS = np.zeros_like(S)
        dotS[:, :self.D] = dotR # Preenche com as velocidades do estado

        # Interação entre os pares de corpos:
        for body_i, body_j in self.pairs:
            '''
            self.pairs: definido no início da simulação
            body_i, body_j: índices dos corpos
            '''

            # Vetor do corpo i até o corpo j
            R_i, R_j = R[body_i], R[body_j] # Posições dos corpos i e j
            Rij = R_j - R_i # Vetor de i => j
            r = np.linalg.norm(Rij) # Módulo do vetor Rij

            # Força de i => j
            F = self.G * self.masses[body_i] * self.masses[body_j] * Rij / (r**2 + self.eps**2)**1.5
            ddotR_i =  F / self.masses[body_i] # Aceleração do corpo i
            ddotR_j = -F / self.masses[body_j] # Aceleração do corpo j

            # Preenche com as acelerações do estado
            dotS[body_i, self.D:] = dotS[body_i, self.D:] + ddotR_i
            dotS[body_j, self.D:] = dotS[body_j, self.D:] + ddotR_j

        return dotS

    def rk4(self, S, dt, func):
        '''
        Entrada:
        - S: Vetor de estado atual
        - dt: Passo de integração
        - func: Função que retorna a derivada do vetor de estado
        Saída:
        - S: Vetor de estado atualizado um passo a frente
        '''

        # Calculando termos do método RK4
        k1 = func(S)
        k2 = func(S + 0.5*k1*dt)
        k3 = func(S + 0.5*k2*dt)
        k4 = func(S + k3*dt)

        # Atualiza S
        S_prime = (1/6.)*(k1 + 2*k2 + 2*k3 + k4)

        return S + S_prime * dt

    def run_simulation(self, t, energy=True, com=True):

        # Função que retorna as posições em relação ao cm:
        def c_m(state, masses):

            cm = np.zeros(self.D)               # Vetor com coordenadas do centro de massa
            r_cm = np.zeros((self.N, self.D))   # Matriz das posições com relação ao centro de massa

            # Nessa etapa eu calculo o vetor do CM:
            for i in range(self.N):
                cm = cm + state[i, :self.D] * masses[i]

            cm = cm / np.sum(masses)
            # fim da etapa

            # Nessa etapa eu calculo cada posição em relação ao CM:
            for i in range(self.N):
                r_cm[i] = state[i, :self.D] - cm
            # fim da etapa

            return r_cm

        # Função que salva os arquivos do estado:
        def save_state(S, t, i, energy, com):
            # Salva o novo estado
            pd.DataFrame(S).to_csv('{}history_{:04d}.csv'.format(self.dir, i), index=False)

            if energy:
                # Atualizando arquivo de energias:
                KE, PE = self.get_energy(S) # Consegue a energia do novo estado
                energies.write("{:.4f},{},{},{}\n".format(t, KE, PE, KE+PE)) # Salva a energia

            if com:
                # Salvando posições no referêncial do CM:
                pd.DataFrame(c_m(S, self.masses)).to_csv('{}positions_{:04d}.csv'.format(self.dir, i), index=False)

        # Número de iterações
        n_iter = len(t)

        # Intervalo entre snapshots
        iterBetSnapshot = int(n_iter*self.nSnp/max(t))

        # Passo de cada iteração
        dt = t[1] - t[0]

        # Determina os índices dos pares de corpos
        self.pairs = list(combinations(range(self.N), 2))

        '''Inicialização:'''

        # Inicializando arquivos de históricos
        # Preenche com as posições iniciais
        pd.DataFrame(self.Si).to_csv('{}history_0000.csv'.format(self.dir), index=False)

        # Inicializando arquivo de energias
        energies = open('{}Energies.csv'.format(self.dir), 'w')
        KE, PE = self.get_energy(self.Si) # Calculando energias iniciais
        energies.write("t,KE,PE,TE\n{},{},{},{}\n".format(0, KE, PE, KE+PE))

        # Inicializando vetor de posições em relação ao centro de massa
        if com:
            # Preenche com as posições iniciais:
            pd.DataFrame(c_m(self.Si, self.masses)).to_csv('{}positions_0000.csv'.format(self.dir), index=False)

        # Simulação
        S = copy.deepcopy(self.Si) # Copia o vetor de estado inicial sem modificá-lo
        info = open('{}info.txt'.format(self.dir), 'w') # Arquivo com informações ao longo da simulação

        startIntegration = time.time() # Tempo inicial do sistema:

        for i in range(1,n_iter+1):
            startStep = time.time()

            # Integração numérica:
            S = self.rk4(S, dt, self.get_state_deriv) # Consegue um novo estado

            # Salvando as snapshots do estado:
            if not i%iterBetSnapshot:
                save_state(S, i*dt, i//iterBetSnapshot, energy, com)

            # Imprimindo as informações do passo de integração:
            if not i%1:
                end = time.time()
                print('Running... TimeSys: {:.8f}, DeltaTimeSys: {:.8f}, TimeSim: {:.4f}, Step: ({}/{})\n'.format(end-startIntegration, end-startStep, i*dt, i, n_iter))
                info.write('Running... TimeSys: {:.8f}, DeltaTimeSys: {:.8f}, TimeSim: {:.4f}, Step: ({}/{})\n'.format(end-startIntegration, end-startStep, i*dt, i, n_iter))

        energies.close()
        info.close()
