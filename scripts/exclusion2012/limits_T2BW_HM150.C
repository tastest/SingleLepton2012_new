float getUpperLimit_T2BW_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 18.7;
     if(seff >= 0.05 && seff < 0.10) ul = 19.0;
     if(seff >= 0.10 && seff < 0.15) ul = 19.4;
     if(seff >= 0.15 && seff < 0.20) ul = 19.7;
     if(seff >= 0.20 && seff < 0.25) ul = 20.3;
     if(seff >= 0.25 && seff < 0.30) ul = 21.0;
     if(seff >= 0.30 && seff < 0.35) ul = 21.8;
     if(seff >= 0.35 && seff < 0.40) ul = 22.5;
     if(seff >= 0.40 && seff < 0.45) ul = 23.4;
     if(seff >= 0.45 && seff < 0.50) ul = 24.4;
     if(seff >= 0.50 && seff < 0.55) ul = 25.3;
     return ul;
}
float getExpectedUpperLimit_T2BW_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 17.9;
     if(seff >= 0.05 && seff < 0.10) ul = 18.0;
     if(seff >= 0.10 && seff < 0.15) ul = 18.4;
     if(seff >= 0.15 && seff < 0.20) ul = 18.6;
     if(seff >= 0.20 && seff < 0.25) ul = 19.2;
     if(seff >= 0.25 && seff < 0.30) ul = 19.9;
     if(seff >= 0.30 && seff < 0.35) ul = 20.5;
     if(seff >= 0.35 && seff < 0.40) ul = 21.3;
     if(seff >= 0.40 && seff < 0.45) ul = 22.1;
     if(seff >= 0.45 && seff < 0.50) ul = 23.0;
     if(seff >= 0.50 && seff < 0.55) ul = 23.9;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 24.7;
     if(seff >= 0.05 && seff < 0.10) ul = 24.8;
     if(seff >= 0.10 && seff < 0.15) ul = 25.4;
     if(seff >= 0.15 && seff < 0.20) ul = 26.4;
     if(seff >= 0.20 && seff < 0.25) ul = 27.6;
     if(seff >= 0.25 && seff < 0.30) ul = 29.0;
     if(seff >= 0.30 && seff < 0.35) ul = 30.7;
     if(seff >= 0.35 && seff < 0.40) ul = 32.6;
     if(seff >= 0.40 && seff < 0.45) ul = 34.7;
     if(seff >= 0.45 && seff < 0.50) ul = 36.8;
     if(seff >= 0.50 && seff < 0.55) ul = 39.3;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 13.2;
     if(seff >= 0.05 && seff < 0.10) ul = 13.2;
     if(seff >= 0.10 && seff < 0.15) ul = 13.3;
     if(seff >= 0.15 && seff < 0.20) ul = 13.6;
     if(seff >= 0.20 && seff < 0.25) ul = 13.9;
     if(seff >= 0.25 && seff < 0.30) ul = 14.1;
     if(seff >= 0.30 && seff < 0.35) ul = 14.3;
     if(seff >= 0.35 && seff < 0.40) ul = 14.7;
     if(seff >= 0.40 && seff < 0.45) ul = 15.0;
     if(seff >= 0.45 && seff < 0.50) ul = 15.5;
     if(seff >= 0.50 && seff < 0.55) ul = 15.9;
     return ul;
}
