<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Pharmacocinétique du baclofène</title>
</head>

<body>
<?php
$TpsMax = 1.5; // Maximum d'assimilatiion en heure
$TpsElim = 3.5; // Demi-vie en hreure
$p = 50;
$yy = 0.50;
$Tmaxp = $TpsMax * 10 * $p; // Tmax assimilation
$Tep = $TpsElim * 10 * $p; // Te émimination

//$ke=SI(Tep;0,5^(1/Tep);1)
if ( $Tep ) {
  $ke = 0.5 ** ( 1 / $Tep );
} else {
  $ke = 1;
}
$yy = $Tmaxp * log( 2 ) / $Tep;
$a = log( 2 ) * $Tmaxp;
$b = log( $Tep ) + ( $a / $Tep );

//CE3 (ta)   =SI(ET(yy>0,99;yy<=1);Tep*0,99;SI(yy>1;Tmaxp*2*LN(2);Tmaxp/10))

if ( $yy > 0.99 and $yy <= 1 ) {
  $ta = $Tep * 0.99;
} elseif ( $yy > 1 ) {
  $ta = $Tmaxp2 * 2 * log( 2 );
} else {
  $ta = $Tmaxp / 10;
}

/*$CO3 = $CE3 / ( $CE3 - $a ) * ( $CE3 * ( $b + 1 - log( $CE3 ) ) - ( 2 * $a ) );
$CY3 = $CO3 / ( $CO3 - $a ) * ( $CO3 * ( $b + 1 - log( $CO3 ) ) - ( 2 * $a ) );
$CY4 = $CY3 / ( $CY3 - $a ) * ( $CY3 * ( $b + 1 - log( $CY3 ) ) - ( 2 * $a ) );
$DI4 = $CY4 / ( $CY4 - $a ) * ( $CY4 * ( $b + 1 - log( $CY4 ) ) - ( 2 * $a ) );
$DI5 = $DI4 / ( $DI4 - $a ) * ( $DI4 * ( $b + 1 - log( $DI4 ) ) - ( 2 * $a ) );
$DS5 = $DI5 / ( $DI5 - $a ) * ( $DI5 * ( $b + 1 - log( $DI5 ) ) - ( 2 * $a ) );
$DS6 = $DS5 / ( $DS5 - $a ) * ( $DS5 * ( $b + 1 - log( $DS5 ) ) - ( 2 * $a ) );
$EC6 = $DS6 / ( $DS6 - $a ) * ( $DS6 * ( $b + 1 - log( $DS6 ) ) - ( 2 * $a ) );
$EC7 = $EC6 / ( $EC6 - $a ) * ( $EC6 * ( $b + 1 - log( $EC6 ) ) - ( 2 * $a ) );
$EM7 = $EC7 / ( $EC7 - $a ) * ( $EC7 * ( $b + 1 - log( $EC7 ) ) - ( 2 * $a ) );
$EM8 = $EM7 / ( $EM7 - $a ) * ( $EM7 * ( $b + 1 - log( $EM7 ) ) - ( 2 * $a ) );
$EW8 = $EM8 / ( $EM8 - $a ) * ( $EM8 * ( $b + 1 - log( $EM8 ) ) - ( 2 * $a ) );
$EW9 = $EW8 / ( $EW8 - $a ) * ( $EW8 * ( $b + 1 - log( $EW8 ) ) - ( 2 * $a ) );
$FG9 = $EW9 / ( $EW9 - $a ) * ( $EW9 * ( $b + 1 - log( $EW9 ) ) - ( 2 * $a ) );
$FG10 = $FG9 / ( $FG9 - $a ) * ( $FG9 * ( $b + 1 - log( $FG9 ) ) - ( 2 * $a ) );
$FG13 = $FG10 / ( $FG10 - $a ) * ( $FG10 * ( $b + 1 - log( $FG10 ) ) - ( 2 * $a ) ); // Ta absorbtion*/

for ( $i = 0; $i < 17; $i++ ) {
  $ta = $ta / ( $ta - $a ) * ( $ta * ( $b + 1 - log( $ta ) ) - ( 2 * $a ) ); // Ta absorbtion
}

$ka0 = 1 / 2 ** ( 1 / $ta );
$ka1 = 0;
//$ka1   =SI(Tmaxp=0;0;SI(OU(Tep=0;Tmaxp<0);0,5^(1/ABS(Tmaxp));ka0+SI(ka0=ke;0,00000000000001)))
if ( $Tmaxp == 0 ) {
  $ka1 = 0;
} elseif ( $Tep == 0 or $Tmaxp < 0 ) {
  $ka1 = 0.5 ** ( 1 / abs( $Tmaxp ) );
} elseif ( $ka0 == $ke ) {
  $ka1 = $ka0 + 0.000000001;
}

// $ka    =SI(Tmaxp=0;0;SI(OU(Tep=0;Tmaxp<0);0,5^(1/ABS(Tmaxp));ka0+SI(ka0=ke;0,00000000000001)))
$ka = $ka1;
if ( $ka1 == $ke ) {
  $ka1 += 0.000000001;
}
// TmamxpVerif 		=SI(OU(ABS(TpsMax)<0,1;ka=1;ke=1);Tmaxp;SI(ke<>ka;LN(-LN(ka)/-LN(ke))/(LN(ke)-LN(ka));Tmaxp))
if ( abs( $Tmaxp ) > 0.1 or $ka == 1 or $ke == 1 ) {
  $TmamxpVerif = $Tmaxp;
} elseif ( $ke <> $ka ) {
  $TmamxpVerif = log( -log( $ka ) / -log( $ke ) ) / ( log( $ke ) - log( $ka ) );
} else {
  $TmamxpVerif = $Tmaxp;
}
echo( "Max d'assimilation : " . $TpsMax . "</br>" );
echo( "Demi-vie d'élimination : " . $TpsElim . "</br>" );
echo( "Te élimination : " . $Tep . "</br>" );
echo( "Ta absorption : " . $ta . "</br>" );
echo( "ke = " . $ke . "</br>" );
echo( "ka = " . $ka . "</br>" );
echo( "Tmax initial : " . $Tmaxp . "</br>" );
echo( "Vérification Tmax : " . $TmamxpVerif . "</br>" );
if ( $Tmaxp = $TmamxpVerif ) {
  echo( "Vérification OK" );
} else {
  echo( "Problèmes de calculs." );
}
?>
</body>
</html>
