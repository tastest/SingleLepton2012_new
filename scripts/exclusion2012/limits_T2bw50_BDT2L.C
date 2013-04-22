float getUpperLimit_T2bw50_BDT2L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.5;
     if(seff >= 0.05 && seff < 0.10) ul = 11.6;
     if(seff >= 0.10 && seff < 0.15) ul = 11.7;
     if(seff >= 0.15 && seff < 0.20) ul = 11.8;
     if(seff >= 0.20 && seff < 0.25) ul = 11.9;
     if(seff >= 0.25 && seff < 0.30) ul = 12.1;
     if(seff >= 0.30 && seff < 0.35) ul = 12.4;
     if(seff >= 0.35 && seff < 0.40) ul = 12.6;
     if(seff >= 0.40 && seff < 0.45) ul = 12.9;
     if(seff >= 0.45 && seff < 0.50) ul = 13.2;
     if(seff >= 0.50 && seff < 0.55) ul = 13.5;
     return ul;
}
float getExpectedUpperLimit_T2bw50_BDT2L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 18.1;
     if(seff >= 0.05 && seff < 0.10) ul = 18.3;
     if(seff >= 0.10 && seff < 0.15) ul = 18.7;
     if(seff >= 0.15 && seff < 0.20) ul = 18.9;
     if(seff >= 0.20 && seff < 0.25) ul = 19.4;
     if(seff >= 0.25 && seff < 0.30) ul = 19.8;
     if(seff >= 0.30 && seff < 0.35) ul = 20.4;
     if(seff >= 0.35 && seff < 0.40) ul = 21.1;
     if(seff >= 0.40 && seff < 0.45) ul = 21.7;
     if(seff >= 0.45 && seff < 0.50) ul = 22.5;
     if(seff >= 0.50 && seff < 0.55) ul = 23.2;
     return ul;
}
float getExpectedP1UpperLimit_T2bw50_BDT2L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 25.1;
     if(seff >= 0.05 && seff < 0.10) ul = 25.5;
     if(seff >= 0.10 && seff < 0.15) ul = 26.1;
     if(seff >= 0.15 && seff < 0.20) ul = 27.0;
     if(seff >= 0.20 && seff < 0.25) ul = 28.3;
     if(seff >= 0.25 && seff < 0.30) ul = 29.7;
     if(seff >= 0.30 && seff < 0.35) ul = 31.2;
     if(seff >= 0.35 && seff < 0.40) ul = 33.0;
     if(seff >= 0.40 && seff < 0.45) ul = 34.9;
     if(seff >= 0.45 && seff < 0.50) ul = 36.9;
     if(seff >= 0.50 && seff < 0.55) ul = 39.0;
     return ul;
}
float getExpectedM1UpperLimit_T2bw50_BDT2L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 13.2;
     if(seff >= 0.05 && seff < 0.10) ul = 13.1;
     if(seff >= 0.10 && seff < 0.15) ul = 13.0;
     if(seff >= 0.15 && seff < 0.20) ul = 13.0;
     if(seff >= 0.20 && seff < 0.25) ul = 13.1;
     if(seff >= 0.25 && seff < 0.30) ul = 13.2;
     if(seff >= 0.30 && seff < 0.35) ul = 13.4;
     if(seff >= 0.35 && seff < 0.40) ul = 13.6;
     if(seff >= 0.40 && seff < 0.45) ul = 13.8;
     if(seff >= 0.45 && seff < 0.50) ul = 13.8;
     if(seff >= 0.50 && seff < 0.55) ul = 14.0;
     return ul;
}
