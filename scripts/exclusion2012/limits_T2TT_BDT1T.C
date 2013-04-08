float getUpperLimit_T2TT_BDT1T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 27.4;
     if(seff >= 0.05 && seff < 0.10) ul = 27.7;
     if(seff >= 0.10 && seff < 0.15) ul = 27.8;
     if(seff >= 0.15 && seff < 0.20) ul = 27.8;
     if(seff >= 0.20 && seff < 0.25) ul = 28.4;
     if(seff >= 0.25 && seff < 0.30) ul = 29.0;
     if(seff >= 0.30 && seff < 0.35) ul = 29.8;
     if(seff >= 0.35 && seff < 0.40) ul = 30.8;
     if(seff >= 0.40 && seff < 0.45) ul = 31.8;
     if(seff >= 0.45 && seff < 0.50) ul = 32.6;
     if(seff >= 0.50 && seff < 0.55) ul = 33.6;
     return ul;
}
float getExpectedUpperLimit_T2TT_BDT1T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 36.0;
     if(seff >= 0.05 && seff < 0.10) ul = 36.6;
     if(seff >= 0.10 && seff < 0.15) ul = 36.9;
     if(seff >= 0.15 && seff < 0.20) ul = 37.9;
     if(seff >= 0.20 && seff < 0.25) ul = 38.8;
     if(seff >= 0.25 && seff < 0.30) ul = 40.0;
     if(seff >= 0.30 && seff < 0.35) ul = 41.3;
     if(seff >= 0.35 && seff < 0.40) ul = 42.7;
     if(seff >= 0.40 && seff < 0.45) ul = 44.2;
     if(seff >= 0.45 && seff < 0.50) ul = 45.8;
     if(seff >= 0.50 && seff < 0.55) ul = 47.5;
     return ul;
}
float getExpectedP1UpperLimit_T2TT_BDT1T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 49.6;
     if(seff >= 0.05 && seff < 0.10) ul = 49.4;
     if(seff >= 0.10 && seff < 0.15) ul = 50.9;
     if(seff >= 0.15 && seff < 0.20) ul = 53.1;
     if(seff >= 0.20 && seff < 0.25) ul = 55.8;
     if(seff >= 0.25 && seff < 0.30) ul = 58.8;
     if(seff >= 0.30 && seff < 0.35) ul = 62.0;
     if(seff >= 0.35 && seff < 0.40) ul = 65.8;
     if(seff >= 0.40 && seff < 0.45) ul = 69.8;
     if(seff >= 0.45 && seff < 0.50) ul = 73.6;
     if(seff >= 0.50 && seff < 0.55) ul = 78.1;
     return ul;
}
float getExpectedM1UpperLimit_T2TT_BDT1T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 26.6;
     if(seff >= 0.05 && seff < 0.10) ul = 26.8;
     if(seff >= 0.10 && seff < 0.15) ul = 27.1;
     if(seff >= 0.15 && seff < 0.20) ul = 27.0;
     if(seff >= 0.20 && seff < 0.25) ul = 27.1;
     if(seff >= 0.25 && seff < 0.30) ul = 27.3;
     if(seff >= 0.30 && seff < 0.35) ul = 27.7;
     if(seff >= 0.35 && seff < 0.40) ul = 28.3;
     if(seff >= 0.40 && seff < 0.45) ul = 28.9;
     if(seff >= 0.45 && seff < 0.50) ul = 29.2;
     if(seff >= 0.50 && seff < 0.55) ul = 30.0;
     return ul;
}
