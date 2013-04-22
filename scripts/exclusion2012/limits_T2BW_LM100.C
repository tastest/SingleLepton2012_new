float getUpperLimit_T2BW_LM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 314.3;
     if(seff >= 0.05 && seff < 0.10) ul = 316.7;
     if(seff >= 0.10 && seff < 0.15) ul = 321.1;
     if(seff >= 0.15 && seff < 0.20) ul = 329.9;
     if(seff >= 0.20 && seff < 0.25) ul = 334.6;
     if(seff >= 0.25 && seff < 0.30) ul = 346.5;
     if(seff >= 0.30 && seff < 0.35) ul = 358.6;
     if(seff >= 0.35 && seff < 0.40) ul = 371.7;
     if(seff >= 0.40 && seff < 0.45) ul = 384.3;
     if(seff >= 0.45 && seff < 0.50) ul = 401.6;
     if(seff >= 0.50 && seff < 0.55) ul = 420.1;
     return ul;
}
float getExpectedUpperLimit_T2BW_LM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 334.0;
     if(seff >= 0.05 && seff < 0.10) ul = 337.8;
     if(seff >= 0.10 && seff < 0.15) ul = 340.1;
     if(seff >= 0.15 && seff < 0.20) ul = 350.8;
     if(seff >= 0.20 && seff < 0.25) ul = 358.1;
     if(seff >= 0.25 && seff < 0.30) ul = 374.3;
     if(seff >= 0.30 && seff < 0.35) ul = 383.7;
     if(seff >= 0.35 && seff < 0.40) ul = 399.2;
     if(seff >= 0.40 && seff < 0.45) ul = 413.5;
     if(seff >= 0.45 && seff < 0.50) ul = 431.6;
     if(seff >= 0.50 && seff < 0.55) ul = 480.1;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_LM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 494.8;
     if(seff >= 0.05 && seff < 0.10) ul = 498.5;
     if(seff >= 0.10 && seff < 0.15) ul = 477.0;
     if(seff >= 0.15 && seff < 0.20) ul = 540.1;
     if(seff >= 0.20 && seff < 0.25) ul = 557.8;
     if(seff >= 0.25 && seff < 0.30) ul = 593.4;
     if(seff >= 0.30 && seff < 0.35) ul = 579.2;
     if(seff >= 0.35 && seff < 0.40) ul = 609.4;
     if(seff >= 0.40 && seff < 0.45) ul = 644.6;
     if(seff >= 0.45 && seff < 0.50) ul = 675.9;
     if(seff >= 0.50 && seff < 0.55) ul = 1862.4;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_LM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 246.2;
     if(seff >= 0.05 && seff < 0.10) ul = 246.9;
     if(seff >= 0.10 && seff < 0.15) ul = 248.1;
     if(seff >= 0.15 && seff < 0.20) ul = 251.1;
     if(seff >= 0.20 && seff < 0.25) ul = 256.3;
     if(seff >= 0.25 && seff < 0.30) ul = 259.7;
     if(seff >= 0.30 && seff < 0.35) ul = 266.2;
     if(seff >= 0.35 && seff < 0.40) ul = 271.6;
     if(seff >= 0.40 && seff < 0.45) ul = 279.4;
     if(seff >= 0.45 && seff < 0.50) ul = 286.8;
     if(seff >= 0.50 && seff < 0.55) ul = 290.1;
     return ul;
}
