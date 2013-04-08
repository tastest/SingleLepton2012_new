float getUpperLimit_T2tt_LM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.5;
     if(seff >= 0.05 && seff < 0.10) ul = 7.6;
     if(seff >= 0.10 && seff < 0.15) ul = 7.8;
     if(seff >= 0.15 && seff < 0.20) ul = 7.7;
     if(seff >= 0.20 && seff < 0.25) ul = 8.0;
     if(seff >= 0.25 && seff < 0.30) ul = 8.1;
     if(seff >= 0.30 && seff < 0.35) ul = 8.2;
     if(seff >= 0.35 && seff < 0.40) ul = 8.5;
     if(seff >= 0.40 && seff < 0.45) ul = 8.7;
     if(seff >= 0.45 && seff < 0.50) ul = 9.0;
     if(seff >= 0.50 && seff < 0.55) ul = 9.5;
     return ul;
}
float getExpectedUpperLimit_T2tt_LM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 8.9;
     if(seff >= 0.05 && seff < 0.10) ul = 9.0;
     if(seff >= 0.10 && seff < 0.15) ul = 9.1;
     if(seff >= 0.15 && seff < 0.20) ul = 9.2;
     if(seff >= 0.20 && seff < 0.25) ul = 9.3;
     if(seff >= 0.25 && seff < 0.30) ul = 9.6;
     if(seff >= 0.30 && seff < 0.35) ul = 10.0;
     if(seff >= 0.35 && seff < 0.40) ul = 10.2;
     if(seff >= 0.40 && seff < 0.45) ul = 10.6;
     if(seff >= 0.45 && seff < 0.50) ul = 11.0;
     if(seff >= 0.50 && seff < 0.55) ul = 11.3;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_LM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 12.3;
     if(seff >= 0.05 && seff < 0.10) ul = 12.5;
     if(seff >= 0.10 && seff < 0.15) ul = 12.7;
     if(seff >= 0.15 && seff < 0.20) ul = 13.5;
     if(seff >= 0.20 && seff < 0.25) ul = 13.8;
     if(seff >= 0.25 && seff < 0.30) ul = 14.4;
     if(seff >= 0.30 && seff < 0.35) ul = 15.3;
     if(seff >= 0.35 && seff < 0.40) ul = 15.6;
     if(seff >= 0.40 && seff < 0.45) ul = 17.0;
     if(seff >= 0.45 && seff < 0.50) ul = 18.1;
     if(seff >= 0.50 && seff < 0.55) ul = 17.6;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_LM300( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.6;
     if(seff >= 0.05 && seff < 0.10) ul = 6.6;
     if(seff >= 0.10 && seff < 0.15) ul = 6.6;
     if(seff >= 0.15 && seff < 0.20) ul = 6.7;
     if(seff >= 0.20 && seff < 0.25) ul = 6.7;
     if(seff >= 0.25 && seff < 0.30) ul = 6.7;
     if(seff >= 0.30 && seff < 0.35) ul = 6.8;
     if(seff >= 0.35 && seff < 0.40) ul = 7.3;
     if(seff >= 0.40 && seff < 0.45) ul = 6.9;
     if(seff >= 0.45 && seff < 0.50) ul = 7.1;
     if(seff >= 0.50 && seff < 0.55) ul = 7.6;
     return ul;
}
