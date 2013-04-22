float getUpperLimit_T2bw50_BDT2T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 9.6;
     if(seff >= 0.05 && seff < 0.10) ul = 9.6;
     if(seff >= 0.10 && seff < 0.15) ul = 9.7;
     if(seff >= 0.15 && seff < 0.20) ul = 9.9;
     if(seff >= 0.20 && seff < 0.25) ul = 10.0;
     if(seff >= 0.25 && seff < 0.30) ul = 10.3;
     if(seff >= 0.30 && seff < 0.35) ul = 10.6;
     if(seff >= 0.35 && seff < 0.40) ul = 11.0;
     if(seff >= 0.40 && seff < 0.45) ul = 11.5;
     if(seff >= 0.45 && seff < 0.50) ul = 12.0;
     if(seff >= 0.50 && seff < 0.55) ul = 12.3;
     return ul;
}
float getExpectedUpperLimit_T2bw50_BDT2T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 9.8;
     if(seff >= 0.05 && seff < 0.10) ul = 9.8;
     if(seff >= 0.10 && seff < 0.15) ul = 9.8;
     if(seff >= 0.15 && seff < 0.20) ul = 10.1;
     if(seff >= 0.20 && seff < 0.25) ul = 10.2;
     if(seff >= 0.25 && seff < 0.30) ul = 10.5;
     if(seff >= 0.30 && seff < 0.35) ul = 10.7;
     if(seff >= 0.35 && seff < 0.40) ul = 11.2;
     if(seff >= 0.40 && seff < 0.45) ul = 11.8;
     if(seff >= 0.45 && seff < 0.50) ul = 12.2;
     if(seff >= 0.50 && seff < 0.55) ul = 12.5;
     return ul;
}
float getExpectedP1UpperLimit_T2bw50_BDT2T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 13.6;
     if(seff >= 0.05 && seff < 0.10) ul = 13.6;
     if(seff >= 0.10 && seff < 0.15) ul = 13.9;
     if(seff >= 0.15 && seff < 0.20) ul = 14.4;
     if(seff >= 0.20 && seff < 0.25) ul = 15.8;
     if(seff >= 0.25 && seff < 0.30) ul = 16.0;
     if(seff >= 0.30 && seff < 0.35) ul = 16.5;
     if(seff >= 0.35 && seff < 0.40) ul = 17.5;
     if(seff >= 0.40 && seff < 0.45) ul = 18.6;
     if(seff >= 0.45 && seff < 0.50) ul = 19.8;
     if(seff >= 0.50 && seff < 0.55) ul = 21.2;
     return ul;
}
float getExpectedM1UpperLimit_T2bw50_BDT2T( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.1;
     if(seff >= 0.05 && seff < 0.10) ul = 7.2;
     if(seff >= 0.10 && seff < 0.15) ul = 7.3;
     if(seff >= 0.15 && seff < 0.20) ul = 7.4;
     if(seff >= 0.20 && seff < 0.25) ul = 7.8;
     if(seff >= 0.25 && seff < 0.30) ul = 7.7;
     if(seff >= 0.30 && seff < 0.35) ul = 7.6;
     if(seff >= 0.35 && seff < 0.40) ul = 7.9;
     if(seff >= 0.40 && seff < 0.45) ul = 8.0;
     if(seff >= 0.45 && seff < 0.50) ul = 8.2;
     if(seff >= 0.50 && seff < 0.55) ul = 8.4;
     return ul;
}
