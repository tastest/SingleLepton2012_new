float getUpperLimit_T2bw50_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 20.1;
     if(seff >= 0.05 && seff < 0.10) ul = 20.3;
     if(seff >= 0.10 && seff < 0.15) ul = 20.7;
     if(seff >= 0.15 && seff < 0.20) ul = 20.6;
     if(seff >= 0.20 && seff < 0.25) ul = 21.1;
     if(seff >= 0.25 && seff < 0.30) ul = 21.7;
     if(seff >= 0.30 && seff < 0.35) ul = 22.2;
     if(seff >= 0.35 && seff < 0.40) ul = 22.9;
     if(seff >= 0.40 && seff < 0.45) ul = 23.7;
     if(seff >= 0.45 && seff < 0.50) ul = 24.3;
     if(seff >= 0.50 && seff < 0.55) ul = 25.4;
     return ul;
}
float getExpectedUpperLimit_T2bw50_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 25.5;
     if(seff >= 0.05 && seff < 0.10) ul = 25.7;
     if(seff >= 0.10 && seff < 0.15) ul = 26.0;
     if(seff >= 0.15 && seff < 0.20) ul = 26.5;
     if(seff >= 0.20 && seff < 0.25) ul = 27.3;
     if(seff >= 0.25 && seff < 0.30) ul = 28.0;
     if(seff >= 0.30 && seff < 0.35) ul = 29.0;
     if(seff >= 0.35 && seff < 0.40) ul = 29.9;
     if(seff >= 0.40 && seff < 0.45) ul = 30.9;
     if(seff >= 0.45 && seff < 0.50) ul = 32.0;
     if(seff >= 0.50 && seff < 0.55) ul = 33.4;
     return ul;
}
float getExpectedP1UpperLimit_T2bw50_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 34.8;
     if(seff >= 0.05 && seff < 0.10) ul = 35.1;
     if(seff >= 0.10 && seff < 0.15) ul = 35.9;
     if(seff >= 0.15 && seff < 0.20) ul = 37.3;
     if(seff >= 0.20 && seff < 0.25) ul = 38.9;
     if(seff >= 0.25 && seff < 0.30) ul = 41.0;
     if(seff >= 0.30 && seff < 0.35) ul = 43.2;
     if(seff >= 0.35 && seff < 0.40) ul = 45.7;
     if(seff >= 0.40 && seff < 0.45) ul = 48.7;
     if(seff >= 0.45 && seff < 0.50) ul = 51.6;
     if(seff >= 0.50 && seff < 0.55) ul = 54.7;
     return ul;
}
float getExpectedM1UpperLimit_T2bw50_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 18.3;
     if(seff >= 0.05 && seff < 0.10) ul = 18.7;
     if(seff >= 0.10 && seff < 0.15) ul = 18.7;
     if(seff >= 0.15 && seff < 0.20) ul = 18.7;
     if(seff >= 0.20 && seff < 0.25) ul = 19.0;
     if(seff >= 0.25 && seff < 0.30) ul = 19.1;
     if(seff >= 0.30 && seff < 0.35) ul = 19.2;
     if(seff >= 0.35 && seff < 0.40) ul = 19.7;
     if(seff >= 0.40 && seff < 0.45) ul = 20.0;
     if(seff >= 0.45 && seff < 0.50) ul = 20.3;
     if(seff >= 0.50 && seff < 0.55) ul = 21.1;
     return ul;
}
