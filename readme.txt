Dati di AMS e Pamela  SOHO e BESS

AMS -> Nuovo Set  [ 28/05/11 - 31/10/2019]
    i flussi di protoni per AMS partona da E circa 500 MeV fino a circa 57 GeV

Pamela -> Ssdc [09/08/06 - 28/01/14]
    i flussi di protoni per PAMELA partona da E circa 105 MeV fino a circa 47 GeV
    
BESS -> date [ 08/2000 - 08/2002  - 12/2004 -  01/2008]

SOHO -> 20 BR [1995 - 2014 ]

Modello fit: N0 fissato a 1

K coeff: standard model  ( kmodel = 1 ) -> K0(t) : K0*Rigidity*Beta;
K coeff: Power law ( kmodel = 2 ) -> K0(t) e a(t)
K coeff: Potgeiter model ( kmodel = 3 ) -> K0(t) e a(t) e b(t)
K coeff: Potgeiter model Rk free ( kmodel = 4 ) -> K0(t) e a(t) e b(t) e R(K)

Opzione per plot delle singole BR con LIS + Flux del modello + dati 
(plot singoli per BR per AMS02 , PAMELA e BESS)

LIS:
  // ---- LIS MODELS ----
  // 0: Corti et al. 2016 ApJ [DEFAULT]
  // 1: Tomassetti et al. 2017 ApJ [OK]
  // 2: Tomassetti et al. 2017 PRD []
  // 3: Tomassetti et al. 2017 ASR
  // 4: REINA




In questa cartella per AMS02 e PAMELA sono salvate le serie temporali di K0 a e b calcolate con
modello Potgeiter usando Rk fisso a 3 GV



PAMELA Potgeiter model valutato con N0 = 0.98 ottenuto da fit con N0 free
NO è fittato per le prime 36 BR (fase descending ciclo 23) in cui resta stabile

Per le correlazioni. Solar Proxy:
SSN smoothed in SSN_convert.root (time format omogeneo)