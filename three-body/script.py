# Biblioteca necessária
import numpy as np

# Função para obter a derivada do vetor de estado:
def get_dotS(S):
    R_1 = S[0, :2]
    R_2 = S[1, :2]
    R_3 = S[2, :2]

    dotR_1 = S[0, 2:]
    dotR_2 = S[1, 2:]
    dotR_3 = S[2, 2:]

    # Constrói a estrutura do Output
    dotS = np.zeros_like(S)

    # Preenche com as velocidades do estado
    dotS[0, :2] = dotR_1
    dotS[1, :2] = dotR_2
    dotS[2, :2] = dotR_3

    # Módulo dos vetores posição
    r12 = np.linalg.norm(R_2 - R_1) # 1 -> 2
    r13 = np.linalg.norm(R_3 - R_1) # 1 -> 3
    r23 = np.linalg.norm(R_3 - R_2) # 2 -> 3

    # Acelerações de cada corpo
    ddotR_1 = G * m_2 * (R_2 - R_1) / r12**3 + G * m_3 * (R_3 - R_1) / r13**3
    ddotR_2 = G * m_1 * (R_1 - R_2) / r12**3 + G * m_3 * (R_3 - R_2) / r23**3
    ddotR_3 = G * m_1 * (R_1 - R_3) / r13**3 + G * m_2 * (R_2 - R_3) / r23**3

    # Preenche com as acelerações do estado
    dotS[0, 2:] = dotS[0, 2:] + ddotR_1
    dotS[1, 2:] = dotS[1, 2:] + ddotR_2
    dotS[2, 2:] = dotS[2, 2:] + ddotR_3

    return dotS

# Função para obter as energias de cada estado:
def get_energy(S):
    R = S[:, :2]        # Submatriz com todas as posições
    dotR = S[:, 2:]     # Submatriz com todas as velocidades
    masses = np.array([m_1, m_2, m_3]).reshape(3,1) # Matriz com as massas

    # Energia Cinética
    KE = 0.5 * np.sum(np.sum(masses * dotR**2 ))

    # Energia Potencial
    PE = 0
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        r = np.linalg.norm(R[i] - R[j])
        PE -= masses[i] * masses[j] / r
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
G, m_1, m_2, m_3 = 1, 1, 1, 1

cte1 = 0.3471128135672417
cte2 = 0.532726851767674
s_0 = np.array([[-1, 0,  cte1  ,   cte2 ],
                [ 1, 0,  cte1  ,   cte2 ],
                [ 0, 0, -2*cte1, -2*cte2]])

# Simulação
t = np.arange(0, 1.5, 0.001)    # Intervalo de simulação
n_iter = len(t)                 # Número de iterações
dt = t[1] - t[0]                # Passo de cada iteração

# Array com as posições de cada corpo (para a plotagem dos gráficos)
positions = np.zeros((n_iter+1, 3, 2))
positions[0] = s_0[:,:2] # Preenche com as posições iniciais

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
    positions[i+1] = S[:,:2]

    # Consegue a energia do novo estado
    KE, PE = get_energy(S)

    # Armazena as energias
    energies[i+1] = np.array([KE, PE, KE+PE], dtype=object)
