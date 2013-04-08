float getUpperLimit_T2tt_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.2;
     if(seff >= 0.05 && seff < 0.10) ul = 4.1;
     if(seff >= 0.10 && seff < 0.15) ul = 4.2;
     if(seff >= 0.15 && seff < 0.20) ul = 4.2;
     if(seff >= 0.20 && seff < 0.25) ul = 4.4;
     if(seff >= 0.25 && seff < 0.30) ul = 4.3;
     if(seff >= 0.30 && seff < 0.35) ul = 4.5;
     if(seff >= 0.35 && seff < 0.40) ul = 4.4;
     if(seff >= 0.40 && seff < 0.45) ul = 4.4;
     if(seff >= 0.45 && seff < 0.50) ul = 4.4;
     if(seff >= 0.50 && seff < 0.55) ul = 4.3;
     return ul;
}
float getExpectedUpperLimit_T2tt_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.1;
     if(seff >= 0.05 && seff < 0.10) ul = 7.2;
     if(seff >= 0.10 && seff < 0.15) ul = 7.4;
     if(seff >= 0.15 && seff < 0.20) ul = 7.5;
     if(seff >= 0.20 && seff < 0.25) ul = 7.6;
     if(seff >= 0.25 && seff < 0.30) ul = 8.0;
     if(seff >= 0.30 && seff < 0.35) ul = 7.9;
     if(seff >= 0.35 && seff < 0.40) ul = 8.0;
     if(seff >= 0.40 && seff < 0.45) ul = 8.2;
     if(seff >= 0.45 && seff < 0.50) ul = 8.6;
     if(seff >= 0.50 && seff < 0.55) ul = 8.7;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.2;
     if(seff >= 0.05 && seff < 0.10) ul = 10.3;
     if(seff >= 0.10 && seff < 0.15) ul = 10.4;
     if(seff >= 0.15 && seff < 0.20) ul = 10.9;
     if(seff >= 0.20 && seff < 0.25) ul = 12.0;
     if(seff >= 0.25 && seff < 0.30) ul = 13.8;
     if(seff >= 0.30 && seff < 0.35) ul = 13.3;
     if(seff >= 0.35 && seff < 0.40) ul = 14.1;
     if(seff >= 0.40 && seff < 0.45) ul = 13.7;
     if(seff >= 0.45 && seff < 0.50) ul = 14.9;
     if(seff >= 0.50 && seff < 0.55) ul = 15.4;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_HM250( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.2;
     if(seff >= 0.05 && seff < 0.10) ul = 5.2;
     if(seff >= 0.10 && seff < 0.15) ul = 5.1;
     if(seff >= 0.15 && seff < 0.20) ul = 5.2;
     if(seff >= 0.20 && seff < 0.25) ul = 5.3;
     if(seff >= 0.25 && seff < 0.30) ul = 5.3;
     if(seff >= 0.30 && seff < 0.35) ul = 5.2;
     if(seff >= 0.35 && seff < 0.40) ul = 5.2;
     if(seff >= 0.40 && seff < 0.45) ul = 5.2;
     if(seff >= 0.45 && seff < 0.50) ul = 5.3;
     if(seff >= 0.50 && seff < 0.55) ul = 5.2;
     return ul;
}
