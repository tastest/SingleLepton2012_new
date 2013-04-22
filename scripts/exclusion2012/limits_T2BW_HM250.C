float getUpperLimit_T2BW_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.0;
     if(seff >= 0.05 && seff < 0.10) ul = 5.0;
     if(seff >= 0.10 && seff < 0.15) ul = 5.0;
     if(seff >= 0.15 && seff < 0.20) ul = 5.0;
     if(seff >= 0.20 && seff < 0.25) ul = 5.1;
     if(seff >= 0.25 && seff < 0.30) ul = 5.3;
     if(seff >= 0.30 && seff < 0.35) ul = 5.4;
     if(seff >= 0.35 && seff < 0.40) ul = 5.6;
     if(seff >= 0.40 && seff < 0.45) ul = 5.5;
     if(seff >= 0.45 && seff < 0.50) ul = 5.7;
     if(seff >= 0.50 && seff < 0.55) ul = 5.9;
     return ul;
}
float getExpectedUpperLimit_T2BW_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.8;
     if(seff >= 0.05 && seff < 0.10) ul = 7.7;
     if(seff >= 0.10 && seff < 0.15) ul = 7.7;
     if(seff >= 0.15 && seff < 0.20) ul = 7.8;
     if(seff >= 0.20 && seff < 0.25) ul = 8.0;
     if(seff >= 0.25 && seff < 0.30) ul = 8.8;
     if(seff >= 0.30 && seff < 0.35) ul = 8.7;
     if(seff >= 0.35 && seff < 0.40) ul = 9.2;
     if(seff >= 0.40 && seff < 0.45) ul = 9.6;
     if(seff >= 0.45 && seff < 0.50) ul = 9.1;
     if(seff >= 0.50 && seff < 0.55) ul = 9.5;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.7;
     if(seff >= 0.05 && seff < 0.10) ul = 11.0;
     if(seff >= 0.10 && seff < 0.15) ul = 11.2;
     if(seff >= 0.15 && seff < 0.20) ul = 11.6;
     if(seff >= 0.20 && seff < 0.25) ul = 12.0;
     if(seff >= 0.25 && seff < 0.30) ul = 13.7;
     if(seff >= 0.30 && seff < 0.35) ul = 13.7;
     if(seff >= 0.35 && seff < 0.40) ul = 14.2;
     if(seff >= 0.40 && seff < 0.45) ul = 14.5;
     if(seff >= 0.45 && seff < 0.50) ul = 15.7;
     if(seff >= 0.50 && seff < 0.55) ul = 16.2;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.7;
     if(seff >= 0.05 && seff < 0.10) ul = 5.4;
     if(seff >= 0.10 && seff < 0.15) ul = 5.4;
     if(seff >= 0.15 && seff < 0.20) ul = 5.5;
     if(seff >= 0.20 && seff < 0.25) ul = 5.5;
     if(seff >= 0.25 && seff < 0.30) ul = 5.8;
     if(seff >= 0.30 && seff < 0.35) ul = 5.7;
     if(seff >= 0.35 && seff < 0.40) ul = 5.8;
     if(seff >= 0.40 && seff < 0.45) ul = 5.7;
     if(seff >= 0.45 && seff < 0.50) ul = 5.8;
     if(seff >= 0.50 && seff < 0.55) ul = 5.9;
     return ul;
}
