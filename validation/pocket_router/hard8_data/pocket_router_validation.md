# PocketRouter Data Validation

Generated: 2026-05-26T20:04:23

## Question

Can sparse residue routing explain hard apo-to-holo pocket adaptation before training a learned router?

## Selected Cases

| test position | original index | metadata RMSD | ligand |
|---:|---:|---:|---|
| 477 | 24500 | 4.0188 | `8pyx__2__1.B__1.J/1.J.sdf` |
| 327 | 24350 | 3.9616 | `8r9u__2__1.B__1.D/1.D.sdf` |
| 390 | 24413 | 3.3143 | `8u6b__1__1.A__1.C/1.C.sdf` |
| 310 | 24333 | 2.9210 | `8ow3__1__1.A__1.B/1.B.sdf` |
| 377 | 24400 | 2.7470 | `8sfu__2__1.B__1.E/1.E.sdf` |
| 342 | 24365 | 2.6896 | `8qn5__3__1.C__1.L/1.L.sdf` |
| 365 | 24388 | 2.6712 | `8pqh__1__1.A__1.B/1.B.sdf` |
| 347 | 24370 | 2.4535 | `8sbv__2__1.B__1.G/1.G.sdf` |

## Aggregate

Mean apo atom RMSD: 1.7027 A

| top-k | router | motion coverage | mobile recall | holo-contact recall | replay RMSD |
|---:|---|---:|---:|---:|---:|
| 4 | distance | 0.103 | 0.111 | 0.406 | 1.5020 |
| 4 | motion_oracle | 0.419 | 0.375 | 0.030 | 1.0276 |
| 4 | contact_oracle | 0.054 | 0.070 | 0.573 | 1.5927 |
| 4 | contact_change_oracle | 0.341 | 0.302 | 0.043 | 1.1666 |
| 4 | random | 0.077 | 0.078 | 0.074 | 1.5704 |
| 8 | distance | 0.156 | 0.184 | 0.704 | 1.4232 |
| 8 | motion_oracle | 0.555 | 0.751 | 0.061 | 0.7868 |
| 8 | contact_oracle | 0.105 | 0.094 | 0.936 | 1.4926 |
| 8 | contact_change_oracle | 0.420 | 0.407 | 0.184 | 1.0196 |
| 8 | random | 0.157 | 0.151 | 0.149 | 1.4369 |
| 12 | distance | 0.212 | 0.255 | 0.897 | 1.3415 |
| 12 | motion_oracle | 0.644 | 1.000 | 0.139 | 0.6426 |
| 12 | contact_oracle | 0.180 | 0.198 | 1.000 | 1.3569 |
| 12 | contact_change_oracle | 0.473 | 0.488 | 0.302 | 0.9225 |
| 12 | random | 0.228 | 0.225 | 0.227 | 1.3115 |
| 16 | distance | 0.263 | 0.292 | 0.952 | 1.2473 |
| 16 | motion_oracle | 0.710 | 1.000 | 0.298 | 0.5252 |
| 16 | contact_oracle | 0.244 | 0.271 | 1.000 | 1.2490 |
| 16 | contact_change_oracle | 0.510 | 0.499 | 0.421 | 0.8612 |
| 16 | random | 0.304 | 0.299 | 0.320 | 1.1855 |

## Interpretation Gate

- Positive route signal: oracle top-k replay RMSD is much lower than apo RMSD, and distance/non-learned routing beats random on motion or contact recall.
- Negative route signal: oracle top-k replay is close to apo RMSD, meaning sparse pocket memory cannot explain the hard motion.
