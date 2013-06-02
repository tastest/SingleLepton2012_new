float getUpperLimit_T2bw50_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 12.7;
     if(seff >= 0.05 && seff < 0.10) ul = 12.8;
     if(seff >= 0.10 && seff < 0.15) ul = 12.9;
     if(seff >= 0.15 && seff < 0.20) ul = 13.1;
     if(seff >= 0.20 && seff < 0.25) ul = 13.3;
     if(seff >= 0.25 && seff < 0.30) ul = 13.8;
     if(seff >= 0.30 && seff < 0.35) ul = 14.0;
     if(seff >= 0.35 && seff < 0.40) ul = 14.4;
     if(seff >= 0.40 && seff < 0.45) ul = 14.8;
     if(seff >= 0.45 && seff < 0.50) ul = 15.2;
     if(seff >= 0.50 && seff < 0.55) ul = 15.8;
     return ul;
}
float getExpectedUpperLimit_T2bw50_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 16.2;
     if(seff >= 0.05 && seff < 0.10) ul = 16.3;
     if(seff >= 0.10 && seff < 0.15) ul = 16.7;
     if(seff >= 0.15 && seff < 0.20) ul = 17.1;
     if(seff >= 0.20 && seff < 0.25) ul = 17.5;
     if(seff >= 0.25 && seff < 0.30) ul = 17.9;
     if(seff >= 0.30 && seff < 0.35) ul = 18.5;
     if(seff >= 0.35 && seff < 0.40) ul = 19.1;
     if(seff >= 0.40 && seff < 0.45) ul = 19.8;
     if(seff >= 0.45 && seff < 0.50) ul = 20.5;
     if(seff >= 0.50 && seff < 0.55) ul = 21.3;
     return ul;
}
float getExpectedP1UpperLimit_T2bw50_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 22.5;
     if(seff >= 0.05 && seff < 0.10) ul = 22.9;
     if(seff >= 0.10 && seff < 0.15) ul = 23.3;
     if(seff >= 0.15 && seff < 0.20) ul = 24.2;
     if(seff >= 0.20 && seff < 0.25) ul = 25.2;
     if(seff >= 0.25 && seff < 0.30) ul = 26.5;
     if(seff >= 0.30 && seff < 0.35) ul = 27.9;
     if(seff >= 0.35 && seff < 0.40) ul = 29.5;
     if(seff >= 0.40 && seff < 0.45) ul = 31.3;
     if(seff >= 0.45 && seff < 0.50) ul = 33.3;
     if(seff >= 0.50 && seff < 0.55) ul = 35.0;
     return ul;
}
float getExpectedM1UpperLimit_T2bw50_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.7;
     if(seff >= 0.05 && seff < 0.10) ul = 11.9;
     if(seff >= 0.10 && seff < 0.15) ul = 11.9;
     if(seff >= 0.15 && seff < 0.20) ul = 12.0;
     if(seff >= 0.20 && seff < 0.25) ul = 12.0;
     if(seff >= 0.25 && seff < 0.30) ul = 12.1;
     if(seff >= 0.30 && seff < 0.35) ul = 12.3;
     if(seff >= 0.35 && seff < 0.40) ul = 12.6;
     if(seff >= 0.40 && seff < 0.45) ul = 12.8;
     if(seff >= 0.45 && seff < 0.50) ul = 13.0;
     if(seff >= 0.50 && seff < 0.55) ul = 13.4;
     return ul;
}
