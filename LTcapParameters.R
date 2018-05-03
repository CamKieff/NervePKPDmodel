#tissue1
consensus_params <- list(KAach = 0.1035, #first-order
KEach = 0.5153,
DVach = 5.7439,
EC50ach = 5.383,
m2max = 0.5,
chemax = 0.5,
IC50m2 = 6,
IC50che = 6,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
init_params <- list(KAach = 3.0207, #inhibition
KEach = 1.2352,
DVach = 5.813,
EC50ach = 5.383,
m2max = 0.9786,
chemax = 0.7053,
IC50m2 = 6.5897,
IC50che = 6.3960,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
freq_list <- c(0.0968, 0.3,0.7,1, 3, 7, 10, 15, 30)
WconDF <- loadNormalizedDF(1, lower = TRUE, dataDF = "cap", normDF = "cap")

#tissue2['
]
consensus_params <- c(KAach = 0.1571, #first-order
KEach = 1.3478,
DVach = 5.9094,
EC50ach = 5.383,
m2max = 0.5,
chemax = 0.5,
IC50m2 = 6,
IC50che = 6,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
init_params <- c(KAach = 1.8256, #inhibition
KEach = 1.7154,
DVach = 5.9126,
EC50ach = 5.383,
m2max = 0.9718,
chemax = 0.3320,
IC50m2 = 6.3091,
IC50che = 8.2514,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
freq_list <- c(0.10343, 0.3,0.7,1, 3, 7, 10, 15, 30)
WconDF <- loadNormalizedDF(2, lower = TRUE, dataDF = "cap", normDF = "cap")

#tissue5
consensus_params <- c(KAach = 0.1050, #first-order
KEach = 0.6771,
DVach = 5.8202,
EC50ach = 5.383,
m2max = 0.5,
chemax = 0.5,
IC50m2 = 6,
IC50che = 6,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
init_params <- c(KAach = 2.2434, #inhibition
KEach = 1.2903,
DVach = 5.9277,
EC50ach = 5.383,
m2max = 0.9786,
chemax = 0.6872,
IC50m2 = 6.4724,
IC50che = 7.4854,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
freq_list <- c(0.10402, 0.3,0.7,1, 3, 7, 10, 15, 30)
WconDF <- loadNormalizedDF(5, lower = TRUE, dataDF = "cap", normDF = "cap")


#tissue7
consensus_params <- c(KAach = 0.239, #first-order
KEach = 0.7195,
DVach = 5.8526,
EC50ach = 5.383,
m2max = 0.9718,
chemax = 0.332,
IC50m2 = 6.3091,
IC50che = 8.2514,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
init_params <- c(KAach = 2.8522, #inhibition
KEach = 0.7651,
DVach = 5.8292,
EC50ach = 5.383,
m2max = 0.9798,
chemax = 0.1475,
IC50m2 = 6.4064,
IC50che = 6.1431,
KAunk = 1, KEunk = 1, DVunk = 7, EC50unk = 5,MAXunk = 1
)
WconDF <- loadNormalizedDF(7, lower = TRUE, dataDF = "cap", normDF = "cap")
