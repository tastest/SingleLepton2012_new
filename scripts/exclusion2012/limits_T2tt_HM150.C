float getUpperLimit_T2tt_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 11.4;
     if(seff >= 0.05 && seff < 0.10) ul = 11.5;
     if(seff >= 0.10 && seff < 0.15) ul = 11.5;
     if(seff >= 0.15 && seff < 0.20) ul = 11.7;
     if(seff >= 0.20 && seff < 0.25) ul = 12.0;
     if(seff >= 0.25 && seff < 0.30) ul = 12.3;
     if(seff >= 0.30 && seff < 0.35) ul = 12.7;
     if(seff >= 0.35 && seff < 0.40) ul = 13.0;
     if(seff >= 0.40 && seff < 0.45) ul = 13.2;
     if(seff >= 0.45 && seff < 0.50) ul = 13.7;
     if(seff >= 0.50 && seff < 0.55) ul = 14.1;
     return ul;
}
float getExpectedUpperLimit_T2tt_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 14.3;
     if(seff >= 0.05 && seff < 0.10) ul = 14.3;
     if(seff >= 0.10 && seff < 0.15) ul = 14.6;
     if(seff >= 0.15 && seff < 0.20) ul = 14.9;
     if(seff >= 0.20 && seff < 0.25) ul = 15.2;
     if(seff >= 0.25 && seff < 0.30) ul = 15.7;
     if(seff >= 0.30 && seff < 0.35) ul = 15.9;
     if(seff >= 0.35 && seff < 0.40) ul = 16.7;
     if(seff >= 0.40 && seff < 0.45) ul = 17.2;
     if(seff >= 0.45 && seff < 0.50) ul = 18.0;
     if(seff >= 0.50 && seff < 0.55) ul = 18.5;
     return ul;
}
float getExpectedP1UpperLimit_T2tt_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 19.7;
     if(seff >= 0.05 && seff < 0.10) ul = 19.9;
     if(seff >= 0.10 && seff < 0.15) ul = 20.4;
     if(seff >= 0.15 && seff < 0.20) ul = 21.0;
     if(seff >= 0.20 && seff < 0.25) ul = 22.0;
     if(seff >= 0.25 && seff < 0.30) ul = 23.0;
     if(seff >= 0.30 && seff < 0.35) ul = 24.8;
     if(seff >= 0.35 && seff < 0.40) ul = 28.5;
     if(seff >= 0.40 && seff < 0.45) ul = 27.4;
     if(seff >= 0.45 && seff < 0.50) ul = 29.0;
     if(seff >= 0.50 && seff < 0.55) ul = 30.7;
     return ul;
}
float getExpectedM1UpperLimit_T2tt_HM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.4;
     if(seff >= 0.05 && seff < 0.10) ul = 10.5;
     if(seff >= 0.10 && seff < 0.15) ul = 10.4;
     if(seff >= 0.15 && seff < 0.20) ul = 10.4;
     if(seff >= 0.20 && seff < 0.25) ul = 10.5;
     if(seff >= 0.25 && seff < 0.30) ul = 10.6;
     if(seff >= 0.30 && seff < 0.35) ul = 11.5;
     if(seff >= 0.35 && seff < 0.40) ul = 11.7;
     if(seff >= 0.40 && seff < 0.45) ul = 11.2;
     if(seff >= 0.45 && seff < 0.50) ul = 11.5;
     if(seff >= 0.50 && seff < 0.55) ul = 11.7;
     return ul;
}
