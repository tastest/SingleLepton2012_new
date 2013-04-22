float getUpperLimit_T2BW_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 17.3;
     if(seff >= 0.05 && seff < 0.10) ul = 17.4;
     if(seff >= 0.10 && seff < 0.15) ul = 17.3;
     if(seff >= 0.15 && seff < 0.20) ul = 17.7;
     if(seff >= 0.20 && seff < 0.25) ul = 18.0;
     if(seff >= 0.25 && seff < 0.30) ul = 18.3;
     if(seff >= 0.30 && seff < 0.35) ul = 18.7;
     if(seff >= 0.35 && seff < 0.40) ul = 19.3;
     if(seff >= 0.40 && seff < 0.45) ul = 19.8;
     if(seff >= 0.45 && seff < 0.50) ul = 20.5;
     if(seff >= 0.50 && seff < 0.55) ul = 21.0;
     return ul;
}
float getExpectedUpperLimit_T2BW_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 23.3;
     if(seff >= 0.05 && seff < 0.10) ul = 23.7;
     if(seff >= 0.10 && seff < 0.15) ul = 24.1;
     if(seff >= 0.15 && seff < 0.20) ul = 24.5;
     if(seff >= 0.20 && seff < 0.25) ul = 25.2;
     if(seff >= 0.25 && seff < 0.30) ul = 25.9;
     if(seff >= 0.30 && seff < 0.35) ul = 26.7;
     if(seff >= 0.35 && seff < 0.40) ul = 27.6;
     if(seff >= 0.40 && seff < 0.45) ul = 28.7;
     if(seff >= 0.45 && seff < 0.50) ul = 29.6;
     if(seff >= 0.50 && seff < 0.55) ul = 30.6;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 32.3;
     if(seff >= 0.05 && seff < 0.10) ul = 32.7;
     if(seff >= 0.10 && seff < 0.15) ul = 33.4;
     if(seff >= 0.15 && seff < 0.20) ul = 34.6;
     if(seff >= 0.20 && seff < 0.25) ul = 36.4;
     if(seff >= 0.25 && seff < 0.30) ul = 38.1;
     if(seff >= 0.30 && seff < 0.35) ul = 40.3;
     if(seff >= 0.35 && seff < 0.40) ul = 42.6;
     if(seff >= 0.40 && seff < 0.45) ul = 45.2;
     if(seff >= 0.45 && seff < 0.50) ul = 47.6;
     if(seff >= 0.50 && seff < 0.55) ul = 50.4;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_LM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 17.1;
     if(seff >= 0.05 && seff < 0.10) ul = 17.3;
     if(seff >= 0.10 && seff < 0.15) ul = 17.2;
     if(seff >= 0.15 && seff < 0.20) ul = 17.3;
     if(seff >= 0.20 && seff < 0.25) ul = 17.4;
     if(seff >= 0.25 && seff < 0.30) ul = 17.6;
     if(seff >= 0.30 && seff < 0.35) ul = 17.8;
     if(seff >= 0.35 && seff < 0.40) ul = 18.2;
     if(seff >= 0.40 && seff < 0.45) ul = 18.5;
     if(seff >= 0.45 && seff < 0.50) ul = 18.8;
     if(seff >= 0.50 && seff < 0.55) ul = 19.2;
     return ul;
}
