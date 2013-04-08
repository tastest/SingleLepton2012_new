float getUpperLimit_T2TT_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.3;
     if(seff >= 0.05 && seff < 0.10) ul = 6.4;
     if(seff >= 0.10 && seff < 0.15) ul = 6.4;
     if(seff >= 0.15 && seff < 0.20) ul = 6.5;
     if(seff >= 0.20 && seff < 0.25) ul = 6.6;
     if(seff >= 0.25 && seff < 0.30) ul = 6.7;
     if(seff >= 0.30 && seff < 0.35) ul = 6.9;
     if(seff >= 0.35 && seff < 0.40) ul = 7.0;
     if(seff >= 0.40 && seff < 0.45) ul = 7.1;
     if(seff >= 0.45 && seff < 0.50) ul = 7.3;
     if(seff >= 0.50 && seff < 0.55) ul = 7.7;
     return ul;
}
float getExpectedUpperLimit_T2TT_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 9.0;
     if(seff >= 0.05 && seff < 0.10) ul = 9.0;
     if(seff >= 0.10 && seff < 0.15) ul = 9.1;
     if(seff >= 0.15 && seff < 0.20) ul = 9.2;
     if(seff >= 0.20 && seff < 0.25) ul = 9.5;
     if(seff >= 0.25 && seff < 0.30) ul = 9.5;
     if(seff >= 0.30 && seff < 0.35) ul = 9.8;
     if(seff >= 0.35 && seff < 0.40) ul = 10.2;
     if(seff >= 0.40 && seff < 0.45) ul = 10.6;
     if(seff >= 0.45 && seff < 0.50) ul = 10.9;
     if(seff >= 0.50 && seff < 0.55) ul = 10.7;
     return ul;
}
float getExpectedP1UpperLimit_T2TT_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 12.5;
     if(seff >= 0.05 && seff < 0.10) ul = 12.5;
     if(seff >= 0.10 && seff < 0.15) ul = 12.9;
     if(seff >= 0.15 && seff < 0.20) ul = 13.2;
     if(seff >= 0.20 && seff < 0.25) ul = 16.2;
     if(seff >= 0.25 && seff < 0.30) ul = 14.4;
     if(seff >= 0.30 && seff < 0.35) ul = 15.3;
     if(seff >= 0.35 && seff < 0.40) ul = 16.9;
     if(seff >= 0.40 && seff < 0.45) ul = 17.2;
     if(seff >= 0.45 && seff < 0.50) ul = 18.2;
     if(seff >= 0.50 && seff < 0.55) ul = 20.0;
     return ul;
}
float getExpectedM1UpperLimit_T2TT_BDT3( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.5;
     if(seff >= 0.05 && seff < 0.10) ul = 6.5;
     if(seff >= 0.10 && seff < 0.15) ul = 6.5;
     if(seff >= 0.15 && seff < 0.20) ul = 6.5;
     if(seff >= 0.20 && seff < 0.25) ul = 6.6;
     if(seff >= 0.25 && seff < 0.30) ul = 6.5;
     if(seff >= 0.30 && seff < 0.35) ul = 6.6;
     if(seff >= 0.35 && seff < 0.40) ul = 6.8;
     if(seff >= 0.40 && seff < 0.45) ul = 6.7;
     if(seff >= 0.45 && seff < 0.50) ul = 6.9;
     if(seff >= 0.50 && seff < 0.55) ul = 7.4;
     return ul;
}
