float getUpperLimit_T2TT_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.2;
     if(seff >= 0.05 && seff < 0.10) ul = 4.3;
     if(seff >= 0.10 && seff < 0.15) ul = 4.4;
     if(seff >= 0.15 && seff < 0.20) ul = 4.4;
     if(seff >= 0.20 && seff < 0.25) ul = 4.5;
     if(seff >= 0.25 && seff < 0.30) ul = 4.5;
     if(seff >= 0.30 && seff < 0.35) ul = 4.7;
     if(seff >= 0.35 && seff < 0.40) ul = 4.7;
     if(seff >= 0.40 && seff < 0.45) ul = 4.9;
     if(seff >= 0.45 && seff < 0.50) ul = 5.1;
     if(seff >= 0.50 && seff < 0.55) ul = 5.2;
     return ul;
}
float getExpectedUpperLimit_T2TT_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.7;
     if(seff >= 0.05 && seff < 0.10) ul = 4.9;
     if(seff >= 0.10 && seff < 0.15) ul = 4.9;
     if(seff >= 0.15 && seff < 0.20) ul = 4.9;
     if(seff >= 0.20 && seff < 0.25) ul = 5.1;
     if(seff >= 0.25 && seff < 0.30) ul = 5.1;
     if(seff >= 0.30 && seff < 0.35) ul = 5.3;
     if(seff >= 0.35 && seff < 0.40) ul = 5.4;
     if(seff >= 0.40 && seff < 0.45) ul = 5.6;
     if(seff >= 0.45 && seff < 0.50) ul = 5.8;
     if(seff >= 0.50 && seff < 0.55) ul = 5.9;
     return ul;
}
float getExpectedP1UpperLimit_T2TT_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.8;
     if(seff >= 0.05 && seff < 0.10) ul = 6.9;
     if(seff >= 0.10 && seff < 0.15) ul = 7.4;
     if(seff >= 0.15 && seff < 0.20) ul = 7.2;
     if(seff >= 0.20 && seff < 0.25) ul = 7.7;
     if(seff >= 0.25 && seff < 0.30) ul = 7.9;
     if(seff >= 0.30 && seff < 0.35) ul = 8.0;
     if(seff >= 0.35 && seff < 0.40) ul = 8.6;
     if(seff >= 0.40 && seff < 0.45) ul = 9.1;
     if(seff >= 0.45 && seff < 0.50) ul = 8.5;
     if(seff >= 0.50 && seff < 0.55) ul = 9.2;
     return ul;
}
float getExpectedM1UpperLimit_T2TT_BDT4( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 3.6;
     if(seff >= 0.05 && seff < 0.10) ul = 3.6;
     if(seff >= 0.10 && seff < 0.15) ul = 3.7;
     if(seff >= 0.15 && seff < 0.20) ul = 3.7;
     if(seff >= 0.20 && seff < 0.25) ul = 3.8;
     if(seff >= 0.25 && seff < 0.30) ul = 3.7;
     if(seff >= 0.30 && seff < 0.35) ul = 3.8;
     if(seff >= 0.35 && seff < 0.40) ul = 3.7;
     if(seff >= 0.40 && seff < 0.45) ul = 3.8;
     if(seff >= 0.45 && seff < 0.50) ul = 4.0;
     if(seff >= 0.50 && seff < 0.55) ul = 4.1;
     return ul;
}
