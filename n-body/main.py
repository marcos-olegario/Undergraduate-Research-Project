# Simulação com arquivos externos de IC

import pandas as pd
import numpy as np

from nbody import nBody

# ------------------- Parâmetros da Simulação ------------------

# Tempo (inicial, final, step):
t = np.arange(0, 4, 0.0009765625)

# Comprimento de amortecimento:
eps = 1

# Constante Gravitacional:
G = 43007.1

# Endereço da pasta:
folderEnd = "Results-simulation/"

# Endereço das condições iniciais:
folderIC = "IC/IC-Gadget.csv"

# Número de Snapshots:
intSnaps = 0.01

# --------------------- Condições Iniciais ---------------------

initialConditions = pd.read_csv(folderIC) # Cuidado com o header!

s_0 = initialConditions.iloc[:, 1:].values
masses = initialConditions.iloc[:,0].values

# ------------------------- Simulação --------------------------

# Criando o objeto:
n_Body = nBody(Si=s_0, masses=masses, eps=eps, G=G, folder=folderEnd, interSnaps = intSnaps)

# Executando a simulação:
n_Body.run_simulation(t,com=False)
