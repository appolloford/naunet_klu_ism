#ifndef __NAUNET_MACROS_H__
#define __NAUNET_MACROS_H__

#include <sunmatrix/sunmatrix_dense.h>

// 
// clang-format off
#define NAUNET_SUCCESS 0
#define NAUNET_FAIL 1

#define MAX_NSYSTEMS 1

#define NELEMENTS 14
#define NSPECIES 669
#define NHEATPROCS 0
#define NCOOLPROCS 0
#define THERMAL (NHEATPROCS || NCOOLPROCS)
#if (NSPECIES + THERMAL)
#define NEQUATIONS (NSPECIES + THERMAL)
#else
#define NEQUATIONS 1
#endif
#define NREACTIONS 8764
// non-zero terms in jacobian matrix, used in sparse matrix
#define NNZ 20813

#define IDX_ELEM_F 0
#define IDX_ELEM_Cl 1
#define IDX_ELEM_P 2
#define IDX_ELEM_Fe 3
#define IDX_ELEM_Mg 4
#define IDX_ELEM_Na 5
#define IDX_ELEM_Si 6
#define IDX_ELEM_S 7
#define IDX_ELEM_N 8
#define IDX_ELEM_O 9
#define IDX_ELEM_He 10
#define IDX_ELEM_C 11
#define IDX_ELEM_GRAIN 12
#define IDX_ELEM_H 13

#define IDX_GFeI 0
#define IDX_GHeI 1
#define IDX_GMgI 2
#define IDX_GNaI 3
#define IDX_GH2S2I 4
#define IDX_GC4H6I 5
#define IDX_GC4SI 6
#define IDX_GCClI 7
#define IDX_GCH3C3NI 8
#define IDX_GCH3C4HI 9
#define IDX_GCH3C5NI 10
#define IDX_GCH3C6HI 11
#define IDX_GCH3C7NI 12
#define IDX_GClOI 13
#define IDX_GHC2PI 14
#define IDX_GHClI 15
#define IDX_GHFI 16
#define IDX_GHNSiI 17
#define IDX_GHPOI 18
#define IDX_GNH2CNI 19
#define IDX_GNO2I 20
#define IDX_GPNI 21
#define IDX_GSiC3HI 22
#define IDX_GSiCH3I 23
#define IDX_GSiO2I 24
#define IDX_GSiSI 25
#define IDX_GFI 26
#define IDX_GH2SiOI 27
#define IDX_GNCCNI 28
#define IDX_GS2I 29
#define IDX_GSiC2HI 30
#define IDX_C2H4CNI 31
#define IDX_GC3SI 32
#define IDX_GC4NI 33
#define IDX_GC4PI 34
#define IDX_GCH2PHI 35
#define IDX_GCH3CHCH2I 36
#define IDX_GCH3COOHI 37
#define IDX_GH2SI 38
#define IDX_GHCNOI 39
#define IDX_GHCPI 40
#define IDX_GHNC3I 41
#define IDX_GHNCOI 42
#define IDX_GHOCNI 43
#define IDX_GHONCI 44
#define IDX_GHS2I 45
#define IDX_GN2OI 46
#define IDX_GPH2I 47
#define IDX_GPOI 48
#define IDX_GSiC2H2I 49
#define IDX_GSiC4I 50
#define IDX_GSiNCI 51
#define IDX_HC2OI 52
#define IDX_GC3PI 53
#define IDX_GC6H6I 54
#define IDX_GH2O2I 55
#define IDX_GHCSiI 56
#define IDX_GPHI 57
#define IDX_GSO2I 58
#define IDX_GSiH4I 59
#define IDX_NCCNCH3II 60
#define IDX_GC2H4CNI 61
#define IDX_GC2SI 62
#define IDX_GC3OI 63
#define IDX_GC4H3I 64
#define IDX_GC7H2I 65
#define IDX_GC8H2I 66
#define IDX_GC9H2I 67
#define IDX_GCH2CNI 68
#define IDX_GCH3CNI 69
#define IDX_GH2CNI 70
#define IDX_GH2CSI 71
#define IDX_GHC2OI 72
#define IDX_GHC5NI 73
#define IDX_GHC7NI 74
#define IDX_GHC9NI 75
#define IDX_GHCCNI 76
#define IDX_GHCSI 77
#define IDX_GSiC3I 78
#define IDX_GSiH3I 79
#define IDX_GSiNI 80
#define IDX_GSiOI 81
#define IDX_HCCNI 82
#define IDX_C3H2OII 83
#define IDX_CFII 84
#define IDX_ClOII 85
#define IDX_GC2H6I 86
#define IDX_GC3H2I 87
#define IDX_GC5H2I 88
#define IDX_GCCPI 89
#define IDX_GCH2CHCCHI 90
#define IDX_GCH3COCH3I 91
#define IDX_GCH3OCH3I 92
#define IDX_GClI 93
#define IDX_GSiCH2I 94
#define IDX_GSiHI 95
#define IDX_GSiH2I 96
#define IDX_H3C3OII 97
#define IDX_HFII 98
#define IDX_SiFII 99
#define IDX_C8H5II 100
#define IDX_C9H5II 101
#define IDX_CH3C3NII 102
#define IDX_GC10H2I 103
#define IDX_GC4H2I 104
#define IDX_GCH2CHCNI 105
#define IDX_GH2CCCI 106
#define IDX_GHC3NI 107
#define IDX_GOCSI 108
#define IDX_H2FII 109
#define IDX_C4H6I 110
#define IDX_C4PII 111
#define IDX_C5NII 112
#define IDX_C6H6II 113
#define IDX_CH2CHCNII 114
#define IDX_COOCH3II 115
#define IDX_FII 116
#define IDX_GC2H5CNI 117
#define IDX_GC2H5OHI 118
#define IDX_GC2OI 119
#define IDX_GC6H2I 120
#define IDX_GC9NI 121
#define IDX_GCO2I 122
#define IDX_GCPI 123
#define IDX_GO2HI 124
#define IDX_GOCNI 125
#define IDX_GSOI 126
#define IDX_GSiCI 127
#define IDX_H2CNOII 128
#define IDX_H2NCOII 129
#define IDX_H2OCNII 130
#define IDX_PC2H4II 131
#define IDX_C7NII 132
#define IDX_C9NII 133
#define IDX_CH3NHII 134
#define IDX_ClOI 135
#define IDX_GC11I 136
#define IDX_GC2NI 137
#define IDX_GCH2CCH2I 138
#define IDX_GCH2NHI 139
#define IDX_GCH3CCHI 140
#define IDX_GCNOI 141
#define IDX_GHNCI 142
#define IDX_GNSI 143
#define IDX_GSiC2I 144
#define IDX_H2C4NII 145
#define IDX_H2CClII 146
#define IDX_H2NOII 147
#define IDX_HClII 148
#define IDX_HNSII 149
#define IDX_HOCII 150
#define IDX_HSOII 151
#define IDX_HSiO2II 152
#define IDX_NH2CNHII 153
#define IDX_C2H5OHII 154
#define IDX_C2OII 155
#define IDX_CH2CHCCHI 156
#define IDX_CH2OHCOII 157
#define IDX_GCH2CCHI 158
#define IDX_GCH2COI 159
#define IDX_GHNOI 160
#define IDX_GN2I 161
#define IDX_HC4NII 162
#define IDX_HCNOHII 163
#define IDX_HNCOHII 164
#define IDX_HOCNII 165
#define IDX_HeHII 166
#define IDX_PC2H3II 167
#define IDX_PNII 168
#define IDX_PNH3II 169
#define IDX_CH3CSII 170
#define IDX_GC5NI 171
#define IDX_GPI 172
#define IDX_H2C7NII 173
#define IDX_H2C9NII 174
#define IDX_H2O2I 175
#define IDX_H3S2II 176
#define IDX_HONCII 177
#define IDX_PCH3II 178
#define IDX_PNH2II 179
#define IDX_C10H3II 180
#define IDX_C2H5CNHII 181
#define IDX_CClII 182
#define IDX_CH2NH2II 183
#define IDX_CH2OHCH2OII 184
#define IDX_CH3C5NHII 185
#define IDX_CH3C7NHII 186
#define IDX_GC2H5I 187
#define IDX_GC3NI 188
#define IDX_GC7NI 189
#define IDX_GCOOCH3I 190
#define IDX_GHCNI 191
#define IDX_H2CNI 192
#define IDX_H2SiOII 193
#define IDX_H3C5NII 194
#define IDX_HCNOII 195
#define IDX_CClI 196
#define IDX_CH3COOHII 197
#define IDX_ClII 198
#define IDX_GNH3I 199
#define IDX_H2ClII 200
#define IDX_H3C9NII 201
#define IDX_H5C2O2II 202
#define IDX_HCOOHII 203
#define IDX_HPNII 204
#define IDX_C2SII 205
#define IDX_C3SII 206
#define IDX_C4NI 207
#define IDX_CH3COOH2II 208
#define IDX_CH3OCH3II 209
#define IDX_GC10I 210
#define IDX_GC10HI 211
#define IDX_GC8HI 212
#define IDX_GCH2OHCOI 213
#define IDX_GCOOHI 214
#define IDX_GO2I 215
#define IDX_H3C7NII 216
#define IDX_H3SiOII 217
#define IDX_HSO2II 218
#define IDX_NH2CNI 219
#define IDX_CH3C3NHII 220
#define IDX_CH3OCH4II 221
#define IDX_GC2H4I 222
#define IDX_GC9I 223
#define IDX_GC9HI 224
#define IDX_GCSI 225
#define IDX_HC7NII 226
#define IDX_HC9NII 227
#define IDX_HNCOII 228
#define IDX_NO2II 229
#define IDX_PC3HII 230
#define IDX_PCH4II 231
#define IDX_SiH4II 232
#define IDX_SiO2I 233
#define IDX_C2H5CNI 234
#define IDX_C2H5OH2II 235
#define IDX_C2H7II 236
#define IDX_CCPII 237
#define IDX_CH2CHCNHII 238
#define IDX_COOCH3I 239
#define IDX_FI 240
#define IDX_GHSI 241
#define IDX_H2POII 242
#define IDX_HCOOCH3II 243
#define IDX_HN2OII 244
#define IDX_N2OII 245
#define IDX_PC2H2II 246
#define IDX_SiC3H2II 247
#define IDX_SiC4II 248
#define IDX_SiC4HII 249
#define IDX_SiH5II 250
#define IDX_SiNCHII 251
#define IDX_C11I 252
#define IDX_C3OII 253
#define IDX_C4H5II 254
#define IDX_C4H7II 255
#define IDX_CH2OHCHOII 256
#define IDX_CH3CHOII 257
#define IDX_GCH2OHCHOI 258
#define IDX_GHCOOHI 259
#define IDX_H2S2I 260
#define IDX_HNSiII 261
#define IDX_HOCSII 262
#define IDX_NSII 263
#define IDX_OCNII 264
#define IDX_PC4HII 265
#define IDX_PH3II 266
#define IDX_SiCH4II 267
#define IDX_C10II 268
#define IDX_C10H2II 269
#define IDX_C2N2II 270
#define IDX_CH3CNII 271
#define IDX_GC2H2I 272
#define IDX_GC5I 273
#define IDX_GC7HI 274
#define IDX_GCH3COI 275
#define IDX_GHCOOCH3I 276
#define IDX_HC4SII 277
#define IDX_HC5NHII 278
#define IDX_SiC2H3II 279
#define IDX_SiC3HI 280
#define IDX_SiCH3II 281
#define IDX_SiNH2II 282
#define IDX_C11II 283
#define IDX_C3NII 284
#define IDX_C4H4II 285
#define IDX_CH2CNII 286
#define IDX_CH3COCH3II 287
#define IDX_CPII 288
#define IDX_GC6I 289
#define IDX_GC8I 290
#define IDX_GCH2OHI 291
#define IDX_GCH3OHI 292
#define IDX_H2S2II 293
#define IDX_HC2PII 294
#define IDX_HC5NII 295
#define IDX_HFI 296
#define IDX_HONCI 297
#define IDX_HSiSII 298
#define IDX_PCH2II 299
#define IDX_SO2II 300
#define IDX_SiNII 301
#define IDX_C9II 302
#define IDX_CH2OHCOI 303
#define IDX_CH3C5NI 304
#define IDX_CH3C7NI 305
#define IDX_GC3HI 306
#define IDX_GC4HI 307
#define IDX_GC5HI 308
#define IDX_GC6HI 309
#define IDX_H2NCII 310
#define IDX_H2SiOI 311
#define IDX_H3CSII 312
#define IDX_HC2OII 313
#define IDX_HCPII 314
#define IDX_SiC4I 315
#define IDX_SiNCI 316
#define IDX_SiNCII 317
#define IDX_CNOI 318
#define IDX_GC2H3I 319
#define IDX_GCH3OI 320
#define IDX_GH2OI 321
#define IDX_HC3SII 322
#define IDX_HCOOH2II 323
#define IDX_HNSiI 324
#define IDX_HPOII 325
#define IDX_HS2II 326
#define IDX_SiC2H2I 327
#define IDX_SiC3II 328
#define IDX_C2NHII 329
#define IDX_C3H6II 330
#define IDX_C6H7II 331
#define IDX_C8H4II 332
#define IDX_CH2OHCHOI 333
#define IDX_ClI 334
#define IDX_GC7I 335
#define IDX_GSiI 336
#define IDX_HC3OII 337
#define IDX_HCNOI 338
#define IDX_SiC3HII 339
#define IDX_SiCH3I 340
#define IDX_C3H7II 341
#define IDX_C4NII 342
#define IDX_C8II 343
#define IDX_C9H4II 344
#define IDX_CH2COII 345
#define IDX_CH2PHI 346
#define IDX_CH3C3NI 347
#define IDX_GC4I 348
#define IDX_GCH3CHOI 349
#define IDX_GH2COI 350
#define IDX_HOCNI 351
#define IDX_HPOI 352
#define IDX_HS2I 353
#define IDX_SiC2II 354
#define IDX_SiC2H2II 355
#define IDX_C6H4II 356
#define IDX_CH3COI 357
#define IDX_CH3COCH4II 358
#define IDX_CH3OHII 359
#define IDX_COOHI 360
#define IDX_GNOI 361
#define IDX_POII 362
#define IDX_C7H4II 363
#define IDX_C7H5II 364
#define IDX_CH3C6HI 365
#define IDX_H2CSII 366
#define IDX_HC9NI 367
#define IDX_HCOOCH3I 368
#define IDX_NCCNHII 369
#define IDX_NO2I 370
#define IDX_PH2I 371
#define IDX_PNI 372
#define IDX_SiC2HI 373
#define IDX_SiC3I 374
#define IDX_C4PI 375
#define IDX_C7II 376
#define IDX_CH2CNI 377
#define IDX_CH3C4HII 378
#define IDX_GCH4I 379
#define IDX_HCPI 380
#define IDX_HClI 381
#define IDX_HNCOI 382
#define IDX_C3PI 383
#define IDX_C9NI 384
#define IDX_CH3C4HI 385
#define IDX_CH3COOHI 386
#define IDX_CH3OCH3I 387
#define IDX_GC2HI 388
#define IDX_GNH2I 389
#define IDX_H2CSI 390
#define IDX_PH2II 391
#define IDX_POI 392
#define IDX_S2I 393
#define IDX_S2II 394
#define IDX_SiCII 395
#define IDX_C3SI 396
#define IDX_C6H5II 397
#define IDX_CH2NHI 398
#define IDX_CH3OI 399
#define IDX_GHCOI 400
#define IDX_HC2PI 401
#define IDX_HCSiI 402
#define IDX_SiCH2I 403
#define IDX_C10HII 404
#define IDX_HCSiII 405
#define IDX_OCSII 406
#define IDX_SiC2I 407
#define IDX_C2H6II 408
#define IDX_C2OI 409
#define IDX_C8HII 410
#define IDX_CH2OHI 411
#define IDX_SiNI 412
#define IDX_C5H5II 413
#define IDX_GNHI 414
#define IDX_GSI 415
#define IDX_HCSI 416
#define IDX_NSI 417
#define IDX_SiC2HII 418
#define IDX_SiH2I 419
#define IDX_C3OI 420
#define IDX_C6HII 421
#define IDX_C9H2I 422
#define IDX_CH2CHCNI 423
#define IDX_CH3CHOHII 424
#define IDX_CH3CNHII 425
#define IDX_CH3COII 426
#define IDX_HC7NI 427
#define IDX_C6II 428
#define IDX_C6H6I 429
#define IDX_C7NI 430
#define IDX_HCOOHI 431
#define IDX_SiH2II 432
#define IDX_SiH3II 433
#define IDX_C10H2I 434
#define IDX_C4II 435
#define IDX_CCPI 436
#define IDX_CPI 437
#define IDX_GC3I 438
#define IDX_O2HI 439
#define IDX_SO2I 440
#define IDX_SiCI 441
#define IDX_C2H5OHI 442
#define IDX_HNOI 443
#define IDX_N2OI 444
#define IDX_SiH3I 445
#define IDX_CH3OH2II 446
#define IDX_C2NII 447
#define IDX_C5II 448
#define IDX_C6H3II 449
#define IDX_C8H3II 450
#define IDX_C9HII 451
#define IDX_C9H3II 452
#define IDX_CH2COI 453
#define IDX_HNC3I 454
#define IDX_OCNI 455
#define IDX_SiCH2II 456
#define IDX_SiHII 457
#define IDX_CH3COCH3I 458
#define IDX_HC3NII 459
#define IDX_C8H2I 460
#define IDX_CO2II 461
#define IDX_GOHI 462
#define IDX_HC5NI 463
#define IDX_SiHI 464
#define IDX_SiH4I 465
#define IDX_C7HII 466
#define IDX_C9H2II 467
#define IDX_PHI 468
#define IDX_C7H2I 469
#define IDX_GCOI 470
#define IDX_HCSII 471
#define IDX_PHII 472
#define IDX_C3H4II 473
#define IDX_C8H2II 474
#define IDX_GCH2I 475
#define IDX_HCO2II 476
#define IDX_C5H3II 477
#define IDX_HC3NHII 478
#define IDX_C5HII 479
#define IDX_CH3CHOI 480
#define IDX_CSII 481
#define IDX_C3H5II 482
#define IDX_C7H3II 483
#define IDX_PII 484
#define IDX_CH3CNI 485
#define IDX_NCCNI 486
#define IDX_C7H2II 487
#define IDX_GC2I 488
#define IDX_C3II 489
#define IDX_H3SII 490
#define IDX_O2HII 491
#define IDX_C6H2II 492
#define IDX_CH4II 493
#define IDX_HNOII 494
#define IDX_H2CCCI 495
#define IDX_C6H2I 496
#define IDX_GCNI 497
#define IDX_C2H5II 498
#define IDX_C5H2I 499
#define IDX_C5H2II 500
#define IDX_CH2CCH2I 501
#define IDX_CNII 502
#define IDX_C2H6I 503
#define IDX_GNI 504
#define IDX_C3H2I 505
#define IDX_C2H5I 506
#define IDX_GOI 507
#define IDX_N2II 508
#define IDX_C4SII 509
#define IDX_HCNII 510
#define IDX_CH3CHCH2I 511
#define IDX_C2II 512
#define IDX_C4HII 513
#define IDX_CNCII 514
#define IDX_HSII 515
#define IDX_NHII 516
#define IDX_O2M 517
#define IDX_CHM 518
#define IDX_H2II 519
#define IDX_OCSI 520
#define IDX_GCH3I 521
#define IDX_NH2II 522
#define IDX_PI 523
#define IDX_SM 524
#define IDX_SiSII 525
#define IDX_C2HM 526
#define IDX_C4HM 527
#define IDX_COII 528
#define IDX_CH3CCHI 529
#define IDX_CH5II 530
#define IDX_HC2SII 531
#define IDX_OHM 532
#define IDX_C10HM 533
#define IDX_C3HM 534
#define IDX_C4SI 535
#define IDX_C6HM 536
#define IDX_C8HM 537
#define IDX_C7HM 538
#define IDX_C9HM 539
#define IDX_C10M 540
#define IDX_C5HM 541
#define IDX_C5NM 542
#define IDX_HSI 543
#define IDX_CM 544
#define IDX_CSI 545
#define IDX_OM 546
#define IDX_C9M 547
#define IDX_CH3OHI 548
#define IDX_HM 549
#define IDX_HC3NI 550
#define IDX_SiSI 551
#define IDX_C3HII 552
#define IDX_C3H3II 553
#define IDX_C8M 554
#define IDX_C2M 555
#define IDX_C2SI 556
#define IDX_C3M 557
#define IDX_C4M 558
#define IDX_C6M 559
#define IDX_H3COII 560
#define IDX_C5M 561
#define IDX_C7M 562
#define IDX_SiOII 563
#define IDX_C2NI 564
#define IDX_GCHI 565
#define IDX_SiOHII 566
#define IDX_NH2I 567
#define IDX_C4H3I 568
#define IDX_C2HII 569
#define IDX_C3H2II 570
#define IDX_C5NI 571
#define IDX_CH2II 572
#define IDX_OHII 573
#define IDX_H2OII 574
#define IDX_C3NM 575
#define IDX_GH2I 576
#define IDX_C3NI 577
#define IDX_CHII 578
#define IDX_CO2I 579
#define IDX_O2II 580
#define IDX_C10HI 581
#define IDX_C9HI 582
#define IDX_SiOI 583
#define IDX_FeII 584
#define IDX_FeI 585
#define IDX_MgII 586
#define IDX_NaII 587
#define IDX_MgI 588
#define IDX_NaI 589
#define IDX_C10I 590
#define IDX_C8HI 591
#define IDX_C7I 592
#define IDX_C9I 593
#define IDX_C8I 594
#define IDX_GCI 595
#define IDX_NHI 596
#define IDX_C5I 597
#define IDX_C6I 598
#define IDX_CNM 599
#define IDX_C4I 600
#define IDX_SOI 601
#define IDX_C7HI 602
#define IDX_C2H4II 603
#define IDX_C4H3II 604
#define IDX_N2HII 605
#define IDX_CH2CCHI 606
#define IDX_SOII 607
#define IDX_CH2CCHII 608
#define IDX_C6HI 609
#define IDX_C5HI 610
#define IDX_H2SII 611
#define IDX_NII 612
#define IDX_C3HI 613
#define IDX_H2COII 614
#define IDX_C3I 615
#define IDX_C4HI 616
#define IDX_NOII 617
#define IDX_OII 618
#define IDX_C4H2I 619
#define IDX_C4H2II 620
#define IDX_HCNHII 621
#define IDX_CH2I 622
#define IDX_HNCI 623
#define IDX_NH3II 624
#define IDX_C2H4I 625
#define IDX_C2H3I 626
#define IDX_NH4II 627
#define IDX_N2I 628
#define IDX_O2I 629
#define IDX_SiII 630
#define IDX_H2SI 631
#define IDX_SII 632
#define IDX_SiI 633
#define IDX_C2H3II 634
#define IDX_NOI 635
#define IDX_HCOI 636
#define IDX_CH4I 637
#define IDX_C2HI 638
#define IDX_H2COI 639
#define IDX_HCNI 640
#define IDX_CHI 641
#define IDX_SI 642
#define IDX_C2I 643
#define IDX_OHI 644
#define IDX_NH3I 645
#define IDX_C2H2II 646
#define IDX_CH3II 647
#define IDX_C2H2I 648
#define IDX_CNI 649
#define IDX_GHI 650
#define IDX_CH3I 651
#define IDX_NI 652
#define IDX_H3OII 653
#define IDX_OI 654
#define IDX_HeII 655
#define IDX_HeI 656
#define IDX_CII 657
#define IDX_H2OI 658
#define IDX_HII 659
#define IDX_CI 660
#define IDX_HCOII 661
#define IDX_H3II 662
#define IDX_COI 663
#define IDX_GRAINM 664
#define IDX_GRAIN0I 665
#define IDX_H2I 666
#define IDX_eM 667
#define IDX_HI 668

#if THERMAL
#define IDX_TGAS NSPECIES
#endif
#define IJth(A, i, j) SM_ELEMENT_D(A, i, j)

#endif