# Baseline de regressão da reconstrução

Referência quantitativa sem corte em sigma (issue #2 do fork
`esilvaju/lorenzetti`), usada para detectar mudanças na reconstrução, como as
do corte em sigma (issues #8–#13).

## Por que um ESD fixo

A digitalização **não é reproduzível**: o ruído eletrônico do `PulseGenerator`
usa `TRandom3 m_rng(0)` (semente aleatória,
`reconstruction/calorimeter/CaloCellBuilder/src/PulseGenerator.cxx:30`), então
duas execuções de `digit_trf.py` sobre o mesmo HIT geram ESDs diferentes.
A reconstrução, por outro lado, é determinística. Por isso o baseline fixa o ESD
(`Zee.ESD.root`, versionado aqui) e só refaz `reco_trf.py` + `ntuple_trf.py`.
Duas reconstruções do mesmo ESD geram tabelas **idênticas bit a bit**.

## Verificar

Com o ambiente carregado (`source build/lzt_setup.sh`):

```bash
python tests/baseline/baseline.py check --outdir /tmp/lzt_baseline
```

Roda a reco e a ntuple sobre `Zee.ESD.root`, extrai as tabelas e compara com as
daqui. Termina com código 0 se tudo bater (~20 s).

Outros subcomandos: `run` (só reco + ntuple), `extract` (tabelas a partir de
AOD + NTUPLE), `compare REF_DIR NEW_DIR`. Use `--help` para ver as opções.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `Zee.ESD.root` | entrada fixa da reco (10 eventos Zee, 1 MB) |
| `Zee.EVT.root` | eventos gerados, só para proveniência |
| `SHA256SUMS` | checksums de EVT, HIT (não versionado) e ESD |
| `eventos.csv` | uma linha por evento |
| `clusters.csv` | uma linha por cluster |
| `baseline.py` | run / extract / compare / check |

O HIT (`Zee.HIT.root`, 276 kB) não está no repo. Ele só seria necessário para
redigitalizar, e redigitalizar não reproduz o ESD.

## Como o ESD foi produzido

Commit `52780a4`, imagem `lorenzetti/lorenzetti:latest` (ROOT 6.34.00,
Geant4 11.2.2, Pythia 8.313, HepMC3 3.03.00, g++ 11.4), 2026-10-16:

```bash
gen_zee.py --nov 10 --seed 1234 -o Zee.EVT.root
simu_trf.py -i Zee.EVT.root -o Zee.HIT.root -nt 10
digit_trf.py -i Zee.HIT.root -o Zee.ESD.root
```

| Etapa | Semente | Reproduzível |
|---|---|---|
| gen | `--seed 1234` (o padrão 0 usa o relógio) | sim |
| simu | Geant4 com semente fixa 512 (`ComponentAccumulator`) | sim, mesmo com `-nt 10` |
| digit | `TRandom3(0)` no `PulseGenerator` | **não** |

Rodar esses comandos de novo gera EVT e HIT com o mesmo conteúdo, mas não o
mesmo ESD. Os checksums não batem nem para EVT e HIT, porque o ROOT grava data
e UUID no arquivo. O `SHA256SUMS` identifica estes arquivos; ele não serve para
validar uma regeneração.

## Tabelas

**`eventos.csv`** (lido do AOD), uma linha por evento, ordenada por `event`:
`n_cells`, `n_truth_cells`, `n_clusters`, `n_truth_clusters`, `n_electrons`,
`n_rings`, `n_seeds`, `sum_cells_e` e `sum_truth_cells_e` (MeV).

**`clusters.csv`** (lido da NTUPLE), uma linha por cluster, ordenada por
(`event`, η), com `cluster` = índice dentro do evento:

- energias: `cl_e`, `cl_et`, `cl_e0..3`, `cl_ehad1..3`, `cl_etot`,
  `cl_e233/237/277`, `cl_emaxs1/2`, `cl_e2tsts1`;
- shower shapes: `cl_reta`, `cl_rphi`, `cl_rhad`, `cl_rhad1`, `cl_eratio`,
  `cl_f0..3`, `cl_weta2`;
- truth: `truth_cl_e/et/eta/phi`;
- elétron: `el_et/eta/phi` e PID `el_tight/medium/loose/vloose` (0/1);
- rings: `cl_rings_Rings`, `cl_ringsL0` e `truth_cl_rings_TruthRings`, cada um
  com vetor completo (separado por espaço), `_n` e `_sum`.

**Fora do baseline:** `secondR`, `lambdaCenter`, `secondLambda`, `fracMax`,
`lateralMom` e `longitudinalMom`. Em clusters não-forward esses membros do
`CaloCluster` não são inicializados e saem com lixo de memória, diferente a cada
execução (issue #19). Voltam ao baseline quando a #19 for corrigida.

### Resumo

| Evento | Células | Células truth | Clusters | Elétrons | Σ E células (GeV) | Σ E truth (GeV) |
|---|---|---|---|---|---|---|
| 0 | 1642 | 975 | 2 | 2 | 114.77 | 31.42 |
| 1 | 664 | 0 | 1 | 1 | 111.05 | 0.00 |
| 2 | 1640 | 952 | 2 | 2 | 257.45 | 46.91 |
| 3 | 90 | 0 | 1 | 1 | 158.95 | 0.00 |
| 4 | 1980 | 956 | 2 | 2 | 146.74 | 49.42 |
| 5 | 1912 | 1912 | 2 | 2 | 82.49 | 82.54 |
| 6 | 2041 | 960 | 2 | 2 | 85.49 | 51.49 |
| 7 | 338 | 0 | 1 | 1 | 215.41 | 0.00 |
| 8 | 1964 | 1964 | 2 | 2 | 112.25 | 117.41 |
| 9 | 1410 | 988 | 2 | 2 | 202.23 | 33.09 |

Total: 17 clusters, 17 elétrons (10 tight). Rings: 100 por cluster central, 18
no cluster forward (evento 3, η = −2.64). Nenhum cluster com |η| > 1.6 tem
cluster truth associado (`truth_cl_e = 0`, rings truth vazios).

## Tolerâncias

| Tipo | Regra |
|---|---|
| contagens, PID, número de rings | igualdade exata |
| floats (energias, shapes, rings, somas) | \|novo − ref\| ≤ 1e-6 + 1e-5·\|ref\| |

Na mesma máquina a diferença é zero. A margem relativa de 1e-5 é para variações
de compilador, flags ou CPU (os valores são `float`, com ~7 dígitos, e somas
sobre milhares de células acumulam arredondamento). A parte absoluta de 1e-6
cobre valores perto de zero, como `cl_rhad`. São tolerâncias **iniciais**: se
outra máquina falhar só por arredondamento, ajuste aqui e no topo do
`baseline.py`.

Referência: uma redigitalização do mesmo HIT muda a Σ E das células em cerca de
0,7% a 5,3% por evento, ordens de grandeza acima da tolerância. O `check` falha como
deveria nesse caso (testado).

## Atualizar o baseline

Quando uma mudança **intencional** alterar a reconstrução (por exemplo o corte em
sigma ligado), gere as tabelas novas e versione-as num commit próprio que
explique a mudança:

```bash
python tests/baseline/baseline.py run --outdir /tmp/novo
python tests/baseline/baseline.py extract --aod /tmp/novo/Zee.AOD.root \
    --ntuple /tmp/novo/Zee.NTUPLE.root --outdir tests/baseline
```

Se a mudança for na digitalização (issues #4 e #5, por exemplo), o ESD também
precisa ser trocado. Regenere-o com os comandos acima, atualize o
`SHA256SUMS` e as tabelas.

## Armadilhas ao ler estes arquivos

- **Ordem dos eventos:** a simulação multithread grava os eventos fora de ordem
  (entrada 0 do ESD = evento 5). Case por `EventNumber` / `eventNumber`, não pelo
  índice.
- **Contadores xAOD:** sem dicionário, `tree.CaloCellContainer_Cells_` devolve
  lixo (por exemplo 557614848) e `GetLeaf("..._").GetValue()` devolve 0. Use o
  comprimento de uma folha-membro: `tree.GetLeaf("CaloCellContainer_Cells.e").GetLen()`.
- **PyROOT:** mantenha o `TFile` vivo enquanto usa a `TTree`; se o arquivo for
  coletado, a árvore aponta para memória liberada e dá segfault.
