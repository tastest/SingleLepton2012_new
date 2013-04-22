float getUpperLimit_T2bw75_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.3;
     if(seff >= 0.05 && seff < 0.10) ul = 7.1;
     if(seff >= 0.10 && seff < 0.15) ul = 7.1;
     if(seff >= 0.15 && seff < 0.20) ul = 7.3;
     if(seff >= 0.20 && seff < 0.25) ul = 7.5;
     if(seff >= 0.25 && seff < 0.30) ul = 7.7;
     if(seff >= 0.30 && seff < 0.35) ul = 7.8;
     if(seff >= 0.35 && seff < 0.40) ul = 8.0;
     if(seff >= 0.40 && seff < 0.45) ul = 8.3;
     if(seff >= 0.45 && seff < 0.50) ul = 8.4;
     if(seff >= 0.50 && seff < 0.55) ul = 8.5;
     return ul;
}
float getExpectedUpperLimit_T2bw75_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.5;
     if(seff >= 0.05 && seff < 0.10) ul = 11.4;
     if(seff >= 0.10 && seff < 0.15) ul = 11.6;
     if(seff >= 0.15 && seff < 0.20) ul = 11.7;
     if(seff >= 0.20 && seff < 0.25) ul = 12.1;
     if(seff >= 0.25 && seff < 0.30) ul = 12.3;
     if(seff >= 0.30 && seff < 0.35) ul = 12.6;
     if(seff >= 0.35 && seff < 0.40) ul = 13.0;
     if(seff >= 0.40 && seff < 0.45) ul = 13.5;
     if(seff >= 0.45 && seff < 0.50) ul = 13.8;
     if(seff >= 0.50 && seff < 0.55) ul = 14.3;
     return ul;
}
float getExpectedP1UpperLimit_T2bw75_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 15.9;
     if(seff >= 0.05 && seff < 0.10) ul = 16.0;
     if(seff >= 0.10 && seff < 0.15) ul = 16.5;
     if(seff >= 0.15 && seff < 0.20) ul = 17.1;
     if(seff >= 0.20 && seff < 0.25) ul = 20.4;
     if(seff >= 0.25 && seff < 0.30) ul = 18.9;
     if(seff >= 0.30 && seff < 0.35) ul = 19.9;
     if(seff >= 0.35 && seff < 0.40) ul = 21.0;
     if(seff >= 0.40 && seff < 0.45) ul = 22.2;
     if(seff >= 0.45 && seff < 0.50) ul = 23.3;
     if(seff >= 0.50 && seff < 0.55) ul = 24.8;
     return ul;
}
float getExpectedM1UpperLimit_T2bw75_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 8.4;
     if(seff >= 0.05 && seff < 0.10) ul = 8.2;
     if(seff >= 0.10 && seff < 0.15) ul = 8.3;
     if(seff >= 0.15 && seff < 0.20) ul = 8.3;
     if(seff >= 0.20 && seff < 0.25) ul = 8.5;
     if(seff >= 0.25 && seff < 0.30) ul = 8.3;
     if(seff >= 0.30 && seff < 0.35) ul = 8.5;
     if(seff >= 0.35 && seff < 0.40) ul = 8.6;
     if(seff >= 0.40 && seff < 0.45) ul = 8.8;
     if(seff >= 0.45 && seff < 0.50) ul = 8.9;
     if(seff >= 0.50 && seff < 0.55) ul = 9.0;
     return ul;
}
