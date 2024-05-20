# Biblioteca necessária
import numpy as np

# Função para obter a derivada do vetor de estado:
def get_dotS(S):
    R_1 = S[0, :3]
    R_2 = S[1, :3]

    dotR_1 = S[0, 3:]
    dotR_2 = S[1, 3:]

    # Constrói a estrutura do Output
    dotS = np.zeros_like(S)

    # Preenche com as velocidades do estado
    dotS[0, :3] = dotR_1
    dotS[1, :3] = dotR_2

    R12 = R_2 - R_1             # Vetor posição de 1 -> 2
    r = np.linalg.norm(R12)     # Módulo do vetor R12

    F = G * m_1 * m_2 * R12 / r**3 # Força de 1 -> 2
    ddotR_1 = F / m_1              # Aceleração do corpo 1
    ddotR_2 = - F / m_2            # Aceleração do corpo 2

    # Preenche com as acelerações do estado
    dotS[0, 3:] = dotS[0, 3:] + ddotR_1
    dotS[1, 3:] = dotS[1, 3:] + ddotR_2

    return dotS

# Função para obter as energias de cada estado:
def get_energy(S):
    R = S[:, :3]        # Submatriz com todas as posições
    dotR = S[:, 3:]     # Submatriz com todas as velocidades
    masses = np.array([m_1, m_2]).reshape(2,1) # Matriz com as massas

    # Energia Cinética
    KE = 0.5 * np.sum(np.sum(masses * dotR**2 ))

    # Energia Potencial
    PE = 0
    r = np.linalg.norm(R[0] - R[1])
    PE -= masses[0] * masses[1] / r
    PE *= G

    return KE, PE

# Função de integração numérica pelo método RK4:
def rk4(S, dt, func):

    # Calculando termos do método RK4
    k1 = func(S)
    k2 = func(S + 0.5*k1*dt)
    k3 = func(S + 0.5*k2*dt)
    k4 = func(S + k3*dt)

    # Atualiza S
    S_prime = (1/6.)*(k1 + 2*k2 + 2*k3 + k4)

    return S + S_prime * dt

# Condições Iniciais
G = 6.67430e-20  # Constante Gravitacional, [km**3/(kg * s**2)]
s_0 = np.array([[   0,    0,    0,   10,   20,   30],
                [3000,    0,    0,    0,   40,    0]]) # Vetor de estado
                                                       # inicial em [km]
                                                       # e [km/s]

m_1 = 1.0e26  # kg
m_2 = 1.0e26  # kg


# Simulação
t = np.arange(0, 520, 0.005)    # Intervalo de simulação
n_iter = len(t)                 # Número de iterações
dt = t[1] - t[0]                # Passo de cada iteração

# Array com as posições de cada corpo (para a plotagem dos gráficos)
positions = np.zeros((n_iter+1, 2, 3))
positions[0] = s_0[:,:3] # Preenche com as posições iniciais

# Array com as energias de cada instante
energies = np.zeros((n_iter+1, 3))
KE, PE = get_energy(s_0)
energies[0] = np.array([KE, PE, KE+PE], dtype=object) # Preenche com as
                                                      # energias iniciais

# Main loop
S = s_0
for i in range(n_iter):

    # Consegue a energia do novo estado
    S = rk4(S, dt, get_dotS)

    # Armazena as posições
    positions[i+1] = S[:,:3]

    # Consegue a energia do novo estado
    KE, PE = get_energy(S)

    # Armazena as energias
    energies[i+1] = np.array([KE, PE, KE+PE], dtype=object)
    
