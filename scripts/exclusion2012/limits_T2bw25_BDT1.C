float getUpperLimit_T2bw25_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.9;
     if(seff >= 0.05 && seff < 0.10) ul = 12.0;
     if(seff >= 0.10 && seff < 0.15) ul = 12.2;
     if(seff >= 0.15 && seff < 0.20) ul = 12.3;
     if(seff >= 0.20 && seff < 0.25) ul = 12.6;
     if(seff >= 0.25 && seff < 0.30) ul = 12.8;
     if(seff >= 0.30 && seff < 0.35) ul = 13.2;
     if(seff >= 0.35 && seff < 0.40) ul = 13.5;
     if(seff >= 0.40 && seff < 0.45) ul = 13.9;
     if(seff >= 0.45 && seff < 0.50) ul = 14.2;
     if(seff >= 0.50 && seff < 0.55) ul = 14.7;
     return ul;
}
float getExpectedUpperLimit_T2bw25_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 15.3;
     if(seff >= 0.05 && seff < 0.10) ul = 15.5;
     if(seff >= 0.10 && seff < 0.15) ul = 15.8;
     if(seff >= 0.15 && seff < 0.20) ul = 16.1;
     if(seff >= 0.20 && seff < 0.25) ul = 16.5;
     if(seff >= 0.25 && seff < 0.30) ul = 16.9;
     if(seff >= 0.30 && seff < 0.35) ul = 17.5;
     if(seff >= 0.35 && seff < 0.40) ul = 18.1;
     if(seff >= 0.40 && seff < 0.45) ul = 18.8;
     if(seff >= 0.45 && seff < 0.50) ul = 19.3;
     if(seff >= 0.50 && seff < 0.55) ul = 20.1;
     return ul;
}
float getExpectedP1UpperLimit_T2bw25_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 21.4;
     if(seff >= 0.05 && seff < 0.10) ul = 21.8;
     if(seff >= 0.10 && seff < 0.15) ul = 22.2;
     if(seff >= 0.15 && seff < 0.20) ul = 23.0;
     if(seff >= 0.20 && seff < 0.25) ul = 24.0;
     if(seff >= 0.25 && seff < 0.30) ul = 25.1;
     if(seff >= 0.30 && seff < 0.35) ul = 26.6;
     if(seff >= 0.35 && seff < 0.40) ul = 28.0;
     if(seff >= 0.40 && seff < 0.45) ul = 29.7;
     if(seff >= 0.45 && seff < 0.50) ul = 31.5;
     if(seff >= 0.50 && seff < 0.55) ul = 33.1;
     return ul;
}
float getExpectedM1UpperLimit_T2bw25_BDT1( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.1;
     if(seff >= 0.05 && seff < 0.10) ul = 11.1;
     if(seff >= 0.10 && seff < 0.15) ul = 11.1;
     if(seff >= 0.15 && seff < 0.20) ul = 11.2;
     if(seff >= 0.20 && seff < 0.25) ul = 11.3;
     if(seff >= 0.25 && seff < 0.30) ul = 11.4;
     if(seff >= 0.30 && seff < 0.35) ul = 11.6;
     if(seff >= 0.35 && seff < 0.40) ul = 11.7;
     if(seff >= 0.40 && seff < 0.45) ul = 12.0;
     if(seff >= 0.45 && seff < 0.50) ul = 12.3;
     if(seff >= 0.50 && seff < 0.55) ul = 12.5;
     return ul;
}
