float getUpperLimit_T2bw50_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 8.2;
     if(seff >= 0.05 && seff < 0.10) ul = 8.2;
     if(seff >= 0.10 && seff < 0.15) ul = 8.3;
     if(seff >= 0.15 && seff < 0.20) ul = 8.4;
     if(seff >= 0.20 && seff < 0.25) ul = 8.6;
     if(seff >= 0.25 && seff < 0.30) ul = 8.8;
     if(seff >= 0.30 && seff < 0.35) ul = 8.9;
     if(seff >= 0.35 && seff < 0.40) ul = 9.1;
     if(seff >= 0.40 && seff < 0.45) ul = 9.4;
     if(seff >= 0.45 && seff < 0.50) ul = 9.7;
     if(seff >= 0.50 && seff < 0.55) ul = 10.0;
     return ul;
}
float getExpectedUpperLimit_T2bw50_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.4;
     if(seff >= 0.05 && seff < 0.10) ul = 10.4;
     if(seff >= 0.10 && seff < 0.15) ul = 10.6;
     if(seff >= 0.15 && seff < 0.20) ul = 10.6;
     if(seff >= 0.20 && seff < 0.25) ul = 11.1;
     if(seff >= 0.25 && seff < 0.30) ul = 11.1;
     if(seff >= 0.30 && seff < 0.35) ul = 11.5;
     if(seff >= 0.35 && seff < 0.40) ul = 11.3;
     if(seff >= 0.40 && seff < 0.45) ul = 12.3;
     if(seff >= 0.45 && seff < 0.50) ul = 12.5;
     if(seff >= 0.50 && seff < 0.55) ul = 13.0;
     return ul;
}
float getExpectedP1UpperLimit_T2bw50_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 14.4;
     if(seff >= 0.05 && seff < 0.10) ul = 14.5;
     if(seff >= 0.10 && seff < 0.15) ul = 16.2;
     if(seff >= 0.15 && seff < 0.20) ul = 17.0;
     if(seff >= 0.20 && seff < 0.25) ul = 17.9;
     if(seff >= 0.25 && seff < 0.30) ul = 18.5;
     if(seff >= 0.30 && seff < 0.35) ul = 17.6;
     if(seff >= 0.35 && seff < 0.40) ul = 18.3;
     if(seff >= 0.40 && seff < 0.45) ul = 19.8;
     if(seff >= 0.45 && seff < 0.50) ul = 20.6;
     if(seff >= 0.50 && seff < 0.55) ul = 22.0;
     return ul;
}
float getExpectedM1UpperLimit_T2bw50_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.5;
     if(seff >= 0.05 && seff < 0.10) ul = 7.4;
     if(seff >= 0.10 && seff < 0.15) ul = 7.7;
     if(seff >= 0.15 && seff < 0.20) ul = 7.7;
     if(seff >= 0.20 && seff < 0.25) ul = 7.8;
     if(seff >= 0.25 && seff < 0.30) ul = 8.0;
     if(seff >= 0.30 && seff < 0.35) ul = 7.9;
     if(seff >= 0.35 && seff < 0.40) ul = 8.0;
     if(seff >= 0.40 && seff < 0.45) ul = 8.0;
     if(seff >= 0.45 && seff < 0.50) ul = 8.3;
     if(seff >= 0.50 && seff < 0.55) ul = 8.5;
     return ul;
}
