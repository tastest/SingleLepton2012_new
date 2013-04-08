float getUpperLimit_T2tt_HM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.0;
     if(seff >= 0.05 && seff < 0.10) ul = 4.0;
     if(seff >= 0.10 && seff < 0.15) ul = 4.0;
     if(seff >= 0.15 && seff < 0.20) ul = 4.0;
     if(seff >= 0.20 && seff < 0.25) ul = 4.1;
     if(seff >= 0.25 && seff < 0.30) ul = 4.1;
     if(seff >= 0.30 && seff < 0.35) ul = 4.2;
     if(seff >= 0.35 && seff < 0.40) ul = 4.3;
     if(seff >= 0.40 && seff < 0.45) ul = 4.4;
     if(seff >= 0.45 && seff < 0.50) ul = 4.5;
     if(seff >= 0.50 && seff < 0.55) ul = 4.6;
     return ul;
}
float getExpectedUpperLimit_T2tt_HM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.6;
     if(seff >= 0.05 && seff < 0.10) ul = 5.6;
     if(seff >= 0.10 && seff < 0.15) ul = 5.5;
     if(seff >= 0.15 && seff < 0.20) ul = 5.8;
     if(seff >= 0.20 && seff < 0.25) ul = 5.6;
     if(seff >= 0.25 && seff < 0.30) ul = 6.2;
     if(seff >= 0.30 && seff < 0.35) ul = 6.1;
     if(seff >= 0.35 && seff < 0.40) ul = 6.4;
     if(seff >= 0.40 && seff < 0.45) ul = 6.6;
     if(seff >= 0.45 && seff < 0.50) ul = 6.6;
     if(seff >= 0.50 && seff < 0.55) ul = 6.9;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_HM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.9;
     if(seff >= 0.05 && seff < 0.10) ul = 7.9;
     if(seff >= 0.10 && seff < 0.15) ul = 8.0;
     if(seff >= 0.15 && seff < 0.20) ul = 8.4;
     if(seff >= 0.20 && seff < 0.25) ul = 8.9;
     if(seff >= 0.25 && seff < 0.30) ul = 9.4;
     if(seff >= 0.30 && seff < 0.35) ul = 9.6;
     if(seff >= 0.35 && seff < 0.40) ul = 10.1;
     if(seff >= 0.40 && seff < 0.45) ul = 10.5;
     if(seff >= 0.45 && seff < 0.50) ul = 10.9;
     if(seff >= 0.50 && seff < 0.55) ul = 11.9;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_HM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.1;
     if(seff >= 0.05 && seff < 0.10) ul = 4.1;
     if(seff >= 0.10 && seff < 0.15) ul = 4.1;
     if(seff >= 0.15 && seff < 0.20) ul = 4.0;
     if(seff >= 0.20 && seff < 0.25) ul = 4.1;
     if(seff >= 0.25 && seff < 0.30) ul = 4.1;
     if(seff >= 0.30 && seff < 0.35) ul = 4.2;
     if(seff >= 0.35 && seff < 0.40) ul = 4.2;
     if(seff >= 0.40 && seff < 0.45) ul = 4.3;
     if(seff >= 0.45 && seff < 0.50) ul = 4.3;
     if(seff >= 0.50 && seff < 0.55) ul = 4.3;
     return ul;
}
