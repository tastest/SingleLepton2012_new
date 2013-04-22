float getUpperLimit_T2BW_HM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 36.2;
     if(seff >= 0.05 && seff < 0.10) ul = 36.4;
     if(seff >= 0.10 && seff < 0.15) ul = 37.0;
     if(seff >= 0.15 && seff < 0.20) ul = 38.3;
     if(seff >= 0.20 && seff < 0.25) ul = 39.5;
     if(seff >= 0.25 && seff < 0.30) ul = 41.2;
     if(seff >= 0.30 && seff < 0.35) ul = 43.0;
     if(seff >= 0.35 && seff < 0.40) ul = 44.8;
     if(seff >= 0.40 && seff < 0.45) ul = 46.8;
     if(seff >= 0.45 && seff < 0.50) ul = 48.8;
     if(seff >= 0.50 && seff < 0.55) ul = 51.0;
     return ul;
}
float getExpectedUpperLimit_T2BW_HM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 28.0;
     if(seff >= 0.05 && seff < 0.10) ul = 28.5;
     if(seff >= 0.10 && seff < 0.15) ul = 29.0;
     if(seff >= 0.15 && seff < 0.20) ul = 29.6;
     if(seff >= 0.20 && seff < 0.25) ul = 30.5;
     if(seff >= 0.25 && seff < 0.30) ul = 31.5;
     if(seff >= 0.30 && seff < 0.35) ul = 32.5;
     if(seff >= 0.35 && seff < 0.40) ul = 33.8;
     if(seff >= 0.40 && seff < 0.45) ul = 35.3;
     if(seff >= 0.45 && seff < 0.50) ul = 36.8;
     if(seff >= 0.50 && seff < 0.55) ul = 38.0;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_HM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 39.0;
     if(seff >= 0.05 && seff < 0.10) ul = 39.1;
     if(seff >= 0.10 && seff < 0.15) ul = 40.1;
     if(seff >= 0.15 && seff < 0.20) ul = 41.9;
     if(seff >= 0.20 && seff < 0.25) ul = 43.6;
     if(seff >= 0.25 && seff < 0.30) ul = 46.1;
     if(seff >= 0.30 && seff < 0.35) ul = 48.8;
     if(seff >= 0.35 && seff < 0.40) ul = 52.0;
     if(seff >= 0.40 && seff < 0.45) ul = 55.2;
     if(seff >= 0.45 && seff < 0.50) ul = 58.3;
     if(seff >= 0.50 && seff < 0.55) ul = 62.3;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_HM100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 20.8;
     if(seff >= 0.05 && seff < 0.10) ul = 21.4;
     if(seff >= 0.10 && seff < 0.15) ul = 21.6;
     if(seff >= 0.15 && seff < 0.20) ul = 22.0;
     if(seff >= 0.20 && seff < 0.25) ul = 22.4;
     if(seff >= 0.25 && seff < 0.30) ul = 23.1;
     if(seff >= 0.30 && seff < 0.35) ul = 23.7;
     if(seff >= 0.35 && seff < 0.40) ul = 24.4;
     if(seff >= 0.40 && seff < 0.45) ul = 25.3;
     if(seff >= 0.45 && seff < 0.50) ul = 25.9;
     if(seff >= 0.50 && seff < 0.55) ul = 26.7;
     return ul;
}
