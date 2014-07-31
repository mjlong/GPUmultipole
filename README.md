multipole
=========
multipole method on GPU with CUDA

What are the branches?
before  -------- Let me call it reference, though it wasn't intended to be

beforesoa ----- Deviate from before only in NeutronInfo AoS to SoA

besoa.ori ------ Successful sort

besoa.ori2 ----- SoA fit[]:tttt,tttt,tttt,....;aaaa,aaaa,aaaa,....;ffff,

besoa.ori3 ----- SoA fit[]: same as ori2; split fitting into 3 loops to fit the data structure, but failed

besoa    ------- SoA fit[]:tttt,aaaa,ffff;tttt,aaaa,ffff;........... latest not fastest version
