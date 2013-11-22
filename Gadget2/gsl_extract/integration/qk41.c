/* integration/qk41.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// #include <config.h>
// #include <gsl/gsl_integration.h>
#include "gsl_integration.h"

/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

static const double xgk[21] =   /* abscissae of the 41-point kronrod rule */
{
  0.998859031588277663838315576545863,
  0.993128599185094924786122388471320,
  0.981507877450250259193342994720217,
  0.963971927277913791267666131197277,
  0.940822633831754753519982722212443,
  0.912234428251325905867752441203298,
  0.878276811252281976077442995113078,
  0.839116971822218823394529061701521,
  0.795041428837551198350638833272788,
  0.746331906460150792614305070355642,
  0.693237656334751384805490711845932,
  0.636053680726515025452836696226286,
  0.575140446819710315342946036586425,
  0.510867001950827098004364050955251,
  0.443593175238725103199992213492640,
  0.373706088715419560672548177024927,
  0.301627868114913004320555356858592,
  0.227785851141645078080496195368575,
  0.152605465240922675505220241022678,
  0.076526521133497333754640409398838,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 20-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 20-point gauss rule */

static const double wg[11] =    /* weights of the 20-point gauss rule */
{
  0.017614007139152118311861962351853,
  0.040601429800386941331039952274932,
  0.062672048334109063569506535187042,
  0.083276741576704748724758143222046,
  0.101930119817240435036750135480350,
  0.118194531961518417312377377711382,
  0.131688638449176626898494499748163,
  0.142096109318382051329298325067165,
  0.149172986472603746787828737001969,
  0.152753387130725850698084331955098
};

static const double wgk[21] =   /* weights of the 41-point kronrod rule */
{
  0.003073583718520531501218293246031,
  0.008600269855642942198661787950102,
  0.014626169256971252983787960308868,
  0.020388373461266523598010231432755,
  0.025882133604951158834505067096153,
  0.031287306777032798958543119323801,
  0.036600169758200798030557240707211,
  0.041668873327973686263788305936895,
  0.046434821867497674720231880926108,
  0.050944573923728691932707670050345,
  0.055195105348285994744832372419777,
  0.059111400880639572374967220648594,
  0.062653237554781168025870122174255,
  0.065834597133618422111563556969398,
  0.068648672928521619345623411885368,
  0.071054423553444068305790361723210,
  0.073030690332786667495189417658913,
  0.074582875400499188986581418362488,
  0.075704497684556674659542775376617,
  0.076377867672080736705502835038061,
  0.076600711917999656445049901530102
};

void
gsl_integration_qk41 (const gsl_function * f, double a, double b,
                      double *result, double *abserr,
                      double *resabs, double *resasc)
{
  double fv1[21], fv2[21];
  gsl_integration_qk (21, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc);
}

