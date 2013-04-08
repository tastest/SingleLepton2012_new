float getUpperLimit_T2TT_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 14.7;
     if(seff >= 0.05 && seff < 0.10) ul = 15.5;
     if(seff >= 0.10 && seff < 0.15) ul = 15.5;
     if(seff >= 0.15 && seff < 0.20) ul = 15.7;
     if(seff >= 0.20 && seff < 0.25) ul = 15.6;
     if(seff >= 0.25 && seff < 0.30) ul = 15.8;
     if(seff >= 0.30 && seff < 0.35) ul = 16.2;
     if(seff >= 0.35 && seff < 0.40) ul = 16.5;
     if(seff >= 0.40 && seff < 0.45) ul = 16.8;
     if(seff >= 0.45 && seff < 0.50) ul = 17.2;
     if(seff >= 0.50 && seff < 0.55) ul = 17.6;
     return ul;
}
float getExpectedUpperLimit_T2TT_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 25.2;
     if(seff >= 0.05 && seff < 0.10) ul = 25.7;
     if(seff >= 0.10 && seff < 0.15) ul = 26.2;
     if(seff >= 0.15 && seff < 0.20) ul = 26.6;
     if(seff >= 0.20 && seff < 0.25) ul = 27.3;
     if(seff >= 0.25 && seff < 0.30) ul = 27.9;
     if(seff >= 0.30 && seff < 0.35) ul = 28.6;
     if(seff >= 0.35 && seff < 0.40) ul = 29.4;
     if(seff >= 0.40 && seff < 0.45) ul = 30.2;
     if(seff >= 0.45 && seff < 0.50) ul = 31.2;
     if(seff >= 0.50 && seff < 0.55) ul = 32.2;
     return ul;
}
float getExpectedP1UpperLimit_T2TT_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 34.7;
     if(seff >= 0.05 && seff < 0.10) ul = 35.2;
     if(seff >= 0.10 && seff < 0.15) ul = 36.0;
     if(seff >= 0.15 && seff < 0.20) ul = 37.5;
     if(seff >= 0.20 && seff < 0.25) ul = 39.2;
     if(seff >= 0.25 && seff < 0.30) ul = 41.2;
     if(seff >= 0.30 && seff < 0.35) ul = 43.6;
     if(seff >= 0.35 && seff < 0.40) ul = 46.2;
     if(seff >= 0.40 && seff < 0.45) ul = 48.9;
     if(seff >= 0.45 && seff < 0.50) ul = 51.6;
     if(seff >= 0.50 && seff < 0.55) ul = 54.8;
     return ul;
}
float getExpectedM1UpperLimit_T2TT_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 18.5;
     if(seff >= 0.05 && seff < 0.10) ul = 18.6;
     if(seff >= 0.10 && seff < 0.15) ul = 18.6;
     if(seff >= 0.15 && seff < 0.20) ul = 18.6;
     if(seff >= 0.20 && seff < 0.25) ul = 18.6;
     if(seff >= 0.25 && seff < 0.30) ul = 18.8;
     if(seff >= 0.30 && seff < 0.35) ul = 18.7;
     if(seff >= 0.35 && seff < 0.40) ul = 19.0;
     if(seff >= 0.40 && seff < 0.45) ul = 19.2;
     if(seff >= 0.45 && seff < 0.50) ul = 19.4;
     if(seff >= 0.50 && seff < 0.55) ul = 19.6;
     return ul;
}
