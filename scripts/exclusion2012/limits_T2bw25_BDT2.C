float getUpperLimit_T2bw25_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.2;
     if(seff >= 0.05 && seff < 0.10) ul = 6.3;
     if(seff >= 0.10 && seff < 0.15) ul = 6.4;
     if(seff >= 0.15 && seff < 0.20) ul = 6.4;
     if(seff >= 0.20 && seff < 0.25) ul = 6.5;
     if(seff >= 0.25 && seff < 0.30) ul = 6.6;
     if(seff >= 0.30 && seff < 0.35) ul = 6.8;
     if(seff >= 0.35 && seff < 0.40) ul = 7.0;
     if(seff >= 0.40 && seff < 0.45) ul = 7.2;
     if(seff >= 0.45 && seff < 0.50) ul = 7.3;
     if(seff >= 0.50 && seff < 0.55) ul = 7.6;
     return ul;
}
float getExpectedUpperLimit_T2bw25_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.7;
     if(seff >= 0.05 && seff < 0.10) ul = 8.0;
     if(seff >= 0.10 && seff < 0.15) ul = 8.1;
     if(seff >= 0.15 && seff < 0.20) ul = 8.1;
     if(seff >= 0.20 && seff < 0.25) ul = 8.2;
     if(seff >= 0.25 && seff < 0.30) ul = 8.4;
     if(seff >= 0.30 && seff < 0.35) ul = 8.7;
     if(seff >= 0.35 && seff < 0.40) ul = 8.9;
     if(seff >= 0.40 && seff < 0.45) ul = 9.2;
     if(seff >= 0.45 && seff < 0.50) ul = 9.5;
     if(seff >= 0.50 && seff < 0.55) ul = 9.9;
     return ul;
}
float getExpectedP1UpperLimit_T2bw25_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.0;
     if(seff >= 0.05 && seff < 0.10) ul = 11.2;
     if(seff >= 0.10 && seff < 0.15) ul = 11.4;
     if(seff >= 0.15 && seff < 0.20) ul = 11.8;
     if(seff >= 0.20 && seff < 0.25) ul = 12.4;
     if(seff >= 0.25 && seff < 0.30) ul = 12.8;
     if(seff >= 0.30 && seff < 0.35) ul = 13.3;
     if(seff >= 0.35 && seff < 0.40) ul = 13.6;
     if(seff >= 0.40 && seff < 0.45) ul = 15.2;
     if(seff >= 0.45 && seff < 0.50) ul = 16.0;
     if(seff >= 0.50 && seff < 0.55) ul = 16.9;
     return ul;
}
float getExpectedM1UpperLimit_T2bw25_BDT2( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.6;
     if(seff >= 0.05 && seff < 0.10) ul = 5.7;
     if(seff >= 0.10 && seff < 0.15) ul = 5.6;
     if(seff >= 0.15 && seff < 0.20) ul = 5.7;
     if(seff >= 0.20 && seff < 0.25) ul = 5.8;
     if(seff >= 0.25 && seff < 0.30) ul = 6.0;
     if(seff >= 0.30 && seff < 0.35) ul = 6.0;
     if(seff >= 0.35 && seff < 0.40) ul = 6.1;
     if(seff >= 0.40 && seff < 0.45) ul = 6.2;
     if(seff >= 0.45 && seff < 0.50) ul = 6.3;
     if(seff >= 0.50 && seff < 0.55) ul = 6.4;
     return ul;
}
