Fit = meu ajuste do modelo dos anéis inclinados para os dados brutos do LITTLE THINGS.
A pasta "Oh et al." contém os dados de velocidade de rotação observada, do gás e do disco para as galáxias tratadas por Oh et al.
A pasta "Nosso ajuste" contém os dados de velocidade de rotação observada obtidas pelo ajuste do modelo do anel inclinado para as mesmas galáxias.
O arquivo "tilted-rings-parameters.csv" contém os parâmetros que utilizei no ajuste, retirados de Oh et al. e do artigo master do LITTLE THINGS.
	AR, DEC, Inclinação e Ângulo de Posição são dados em graus; Vsys (Systemic Velocity ) e dv (Channel Separation) são dados em km/s; D e rmax são Mpc e kpc, respectivamente.

Abaixo, segue uma descrição do conjunto de dados:

Sem informações de estrelas, apenas do gás (Removi estas galáxias das análises):
	DDO43
	DDO47
	
Sem fit:
	F564-V3
	IC1613
	UGC8508
	
Fit ajustou aos dados:
	DDO52
	DDO126
	DDO154
	DDO216

Fit ajustou parcialmente aos dados:
	CVnidwa - até 1,5 kpc
	DDO46 - até 3,0 kpc
	DDO70 - Mesmo formato, "shiftado"
	DDO87 - até 5,5 kpc
	DDO101 - Mesmo formato, "shiftado"
	DDO133 - até 2,5 kpc
	DDO168 - até 3,5 kpc
	DDO210 - até 0,25 kpc
	Haro29 - até 3,8 kpc
	Haro36 - Mesmo formato, "shiftado"
	NGC2366 - até 4,0 kpc
	WLM - até 3,0 kpc

Não ajustou aos dados:
	DDO50
	DDO53
	IC10
	NGC1569
	NGC3738

Em resumo:
Galáxias disponíveis para uso (Somente com os dados de Oh et. al.) (24 galáxias):
	CVnidwa
	DDO 46
	DDO 50
	DDO 52
	DDO 53
	DDO 70
	DDO 87
	DDO 101
	DDO 126
	DDO 133
	DDO 154 
	DDO 168 
	DDO 210 
	DDO 216 
	F564-V3
	Haro 29
	Haro 36 
	IC 10
	IC 1613
	NGC 1569
	NGC 2366
	NGC 3738
	UGC 8508
	WLM 

Galáxias disponíveis para uso (dados de Oh et. al. + meu ajuste) (19 galáxias):
(Removendo as galáxias que o fit não ajustou nem ao menos parcialmente)
	CVnidwa
	DDO 46
	DDO 52
	DDO 70
	DDO 87
	DDO 101
	DDO 126
	DDO 133
	DDO 154 
	DDO 168 
	DDO 210 
	DDO 216 
	F564-V3
	Haro 29
	Haro 36 
	IC 1613
	NGC 2366
	UGC 8508
	WLM 
