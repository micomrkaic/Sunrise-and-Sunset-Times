/*
 * This file is part of Mico's toy RPN Calculator
 *
 * Sun Times calculator is free software: 
 * you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mico's toy RPN Calculator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mico's toy RPN Calculator. If not, see <https://www.gnu.org/licenses/>.
 */

// sun_times.c  (C17, no GCC extensions; with input validation)
// Build: gcc sun_times_v2.c -o sun_times -lm


/*
  | Town                        | Latitude                               | Longitude                               |
  | --------------------------- | -------------------------------------- | --------------------------------------- |
  | Kranj, Slovenia             | ~ 46.2389° N ([latitude.to][1])        | ~ 14.3556° E ([Geodatos][2])            |
  | Duluth, Georgia, USA        | ~ 34.0029° N ([latitude.to][3])        | ~ –84.1446° W ([latitude.to][3])        |
  | Alexandria, Virginia, USA   | ~ 38.8048° N ([Coordinates Finder][4]) | ~ –77.0469° W ([latitude.to][5])        |
  | Durham, North Carolina, USA | ~ 35.9940° N ([latitude.to][6])        | ~ –78.8986° W ([Coordinates Finder][7]) |

  [1]: https://latitude.to/map/si/slovenia/cities/kranj?utm_source=chatgpt.com "GPS coordinates of Kranj, Slovenia. Latitude: 46.2389 Longitude"
  [2]: https://www.geodatos.net/en/coordinates/slovenia/kranj?utm_source=chatgpt.com "Kranj Geographic coordinates - Latitude & longitude - Geodatos"
  [3]: https://latitude.to/map/us/united-states/cities/duluth-georgia?utm_source=chatgpt.com "GPS coordinates of Duluth, Georgia, United States. Latitude"
  [4]: https://www.coordinatesfinder.com/coordinates/5333-alexandria-virginia?utm_source=chatgpt.com "GPS coordinates for Alexandria Virginia | CoordinatesFinder.com"
  [5]: https://latitude.to/map/us/united-states/cities/alexandria-virginia?utm_source=chatgpt.com "GPS coordinates of Alexandria, Virginia, United States. Latitude"
  [6]: https://latitude.to/map/us/united-states/cities/durham?utm_source=chatgpt.com "GPS coordinates of Durham, United States. Latitude: 35.9940 ..."
  [7]: https://www.coordinatesfinder.com/coordinates/335771-durham-nc?utm_source=chatgpt.com "GPS coordinates for Durham NC | CoordinatesFinder.com"
*/

/* TODO as of November 1, 2025

   1. If empty options use the IP determined locality and today's date
   2. Add the day of the week in the output
   3. Make the options for the output: only UTC, only local, both, no header or with header

*/

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ======= math helpers =======
static double deg2rad(double d){ return d * (M_PI/180.0); }
static double rad2deg(double r){ return r * (180.0/M_PI); }
static double norm360(double x){ double y = fmod(x, 360.0); return (y < 0) ? y + 360.0 : y; }
static double norm24(double h){ double y = fmod(h, 24.0); return (y < 0) ? y + 24.0 : y; }
static double clamp(double x, double a, double b){ return x < a ? a : (x > b ? b : x); }

// ======= date → JD =======
// Julian Day at 00:00 UTC for Gregorian date (Fliegel–Van Flandern)
static double jd_utc_midnight(int Y,int M,int D){
  int a = (14 - M)/12;
  int y = Y + 4800 - a;
  int m = M + 12*a - 3;
  return D + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
}

// ======= validation =======
static int is_leap_year(int Y){
  return ((Y % 4 == 0) && (Y % 100 != 0)) || (Y % 400 == 0);
}
static int days_in_month(int Y, int M){
  static const int dm[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  if (M == 2) return 28 + is_leap_year(Y);
  return dm[(M-1 >= 0 && M-1 < 12) ? M-1 : 0];
}
static int validate_inputs(int Y,int M,int D,
                           double lat_deg,double lon_deg,
                           double tz_offset_hours,int refine_iters,
                           char *err, size_t errsz)
{
  if (M < 1 || M > 12) {
    if (err && errsz) snprintf(err, errsz, "Invalid month %d (must be 1..12).", M);
    return 1;
  }
  {
    int dim = days_in_month(Y, M);
    if (D < 1 || D > dim) {
      if (err && errsz) snprintf(err, errsz, "Invalid day %d for %04d-%02d (month has %d days).", D, Y, M, dim);
      return 2;
    }
  }
  if (!(lat_deg >= -90.0 && lat_deg <= 90.0)) {
    if (err && errsz) snprintf(err, errsz, "Latitude %.6f out of range [-90, 90].", lat_deg);
    return 3;
  }
  if (!(lon_deg >= -180.0 && lon_deg <= 180.0)) {
    if (err && errsz) snprintf(err, errsz, "Longitude %.6f out of range [-180, 180].", lon_deg);
    return 4;
  }
  if (!isfinite(tz_offset_hours)) {
    if (err && errsz) snprintf(err, errsz, "Time-zone offset must be finite.");
    return 5;
  }
  if (refine_iters < 0 || refine_iters > 10) {
    if (err && errsz) snprintf(err, errsz, "refine_iters %d out of reasonable range (0..10).", refine_iters);
    return 6;
  }
  if (err && errsz) err[0] = '\0';
  return 0;
}

// ======= solar model (NOAA-style approximations) =======
static double T_centuries(double JD){ return (JD - 2451545.0)/36525.0; } // since J2000.0

static double solar_M(double T){ // mean anomaly (deg)
  return norm360(357.52911 + T*(35999.05029 - 0.0001537*T));
}

static double solar_C(double M_deg, double T){ // equation of center (deg)
  double M = deg2rad(M_deg);
  return (1.914602 - T*(0.004817 + 0.000014*T))*sin(M)
    + (0.019993 - 0.000101*T)*sin(2*M)
    + 0.000289*sin(3*M);
}

static double mean_longitude(double T){ // L0 (deg)
  return norm360(280.46646 + T*(36000.76983 + 0.0003032*T));
}

static double obliquity_deg(double T){
  return 23.439291 - T*(0.013004167 + T*(1.63889e-7 - 5.03611e-7*T));
}

static double solar_lambda(double L0_deg, double C_deg){ // ecliptic longitude (deg)
  return norm360(L0_deg + C_deg);
}

static double solar_decl_deg(double lambda_deg, double eps_deg){
  double s = sin(deg2rad(eps_deg)) * sin(deg2rad(lambda_deg));
  return rad2deg(asin(s));
}

static double equation_of_time_min(double T){
  // NOAA approximation, returns minutes
  double eps = deg2rad(obliquity_deg(T));
  double L0  = deg2rad(mean_longitude(T));
  double e   = 0.016708634 - T*(0.000042037 + 0.0000001267*T);
  double M   = deg2rad(solar_M(T));
  double y   = tan(eps/2.0); y *= y;

  double EoT = y*sin(2*L0) - 2*e*sin(M) + 4*e*y*sin(M)*cos(2*L0)
    - 0.5*y*y*sin(4*L0) - 1.25*e*e*sin(2*M);
  return rad2deg(EoT) * 4.0;
}

// Hour angle (deg) for given solar altitude (deg); NaN if impossible.
static double hour_angle_deg(double lat_deg, double dec_deg, double alt_deg){
  double lat = deg2rad(lat_deg);
  double dec = deg2rad(dec_deg);
  double alt = deg2rad(alt_deg);
  double cosH = (sin(alt) - sin(lat)*sin(dec)) / (cos(lat)*cos(dec));
  cosH = clamp(cosH, -1.0, 1.0); // allow tiny numeric overshoots
  if (cosH < -1.0 + 1e-12 || cosH > 1.0 - 1e-12) {
    if (cosH > 1.0 || cosH < -1.0) return NAN; // truly no event
  }
  return rad2deg(acos(cosH)); // [0..180]
}

// ======= public API =======
typedef struct {
  // Flags
  int have_civil_dawn, have_sunrise, have_sunset, have_civil_dusk;
  int polar_day, polar_night; // best-effort flags

  // UTC times (fractional hours [0,24))
  double civil_dawn_utc, sunrise_utc, sunset_utc, civil_dusk_utc;

  // Local times (UTC + tz_offset_hours, fractional hours [0,24))
  double civil_dawn_local, sunrise_local, sunset_local, civil_dusk_local;
} SunTimes;

// Internal solar quantities at a given JD
static void solar_at_JD(double JD, double *dec_deg_out, double *EoT_min_out){
  double T   = T_centuries(JD);
  double M   = solar_M(T);
  double C   = solar_C(M, T);
  double L0  = mean_longitude(T);
  double lam = solar_lambda(L0, C);
  double eps = obliquity_deg(T);
  if (dec_deg_out)  *dec_deg_out  = solar_decl_deg(lam, eps);
  if (EoT_min_out)  *EoT_min_out  = equation_of_time_min(T);
}

// Compute UTC time from solar noon and hour angle/sign (-1 = rise/dawn, +1 = set/dusk)
static inline double make_time(double solar_noon_utc, double Hdeg, double sign) {
  if (isnan(Hdeg)) return NAN;
  return norm24(solar_noon_utc + sign * (Hdeg / 15.0));
}

// Compute events for date Y-M-D at lat, lon (deg), with fixed tz offset (hours).
// refine_iters: 0 (fast) … 2 (very accurate).
// errbuf (optional) can be NULL; when provided, it receives validation errors.
// Returns 0 on success; nonzero if validation failed.
int compute_sun_times(int Y,int M,int D, double lat_deg, double lon_deg,
                      double tz_offset_hours, int refine_iters, SunTimes *out,
                      char *errbuf, size_t errbuf_sz)
{
  if (!out) return 98;
  if (errbuf && errbuf_sz) errbuf[0] = '\0';

  // Validate inputs
  {
    char local_err[160];
    int v = validate_inputs(Y,M,D,lat_deg,lon_deg,tz_offset_hours,refine_iters,
			    local_err, sizeof local_err);
    if (v != 0) {
      if (errbuf && errbuf_sz) snprintf(errbuf, errbuf_sz, "%s", local_err);
      return v;
    }
  }

  // Clamp refine_iters to 0..2
  if (refine_iters < 0) refine_iters = 0;
  if (refine_iters > 2) refine_iters = 2;

  memset(out, 0, sizeof(*out));

  // 1) Base solar quantities at 00:00 UTC of the calendar date
  double JD0 = jd_utc_midnight(Y,M,D);

  double dec0, EoT0_min;
  solar_at_JD(JD0, &dec0, &EoT0_min);

  // 2) Local solar noon in UTC fractional hours
  double solar_noon_utc = 12.0 + (-lon_deg)*4.0/60.0 - EoT0_min/60.0;

  // 3) Hour angles for civil (−6°) and sunrise/set (−0.833°)
  double H_civil_deg = hour_angle_deg(lat_deg, dec0, -6.0);
  double H_sun_deg   = hour_angle_deg(lat_deg, dec0, -0.833);

  // Initial UTC estimates
  double civil_dawn_utc = make_time(solar_noon_utc, H_civil_deg, -1.0);
  double civil_dusk_utc = make_time(solar_noon_utc, H_civil_deg, +1.0);
  double sunrise_utc    = make_time(solar_noon_utc, H_sun_deg,   -1.0);
  double sunset_utc     = make_time(solar_noon_utc, H_sun_deg,   +1.0);

  // 4) Optional refinement iterations
  for (int it = 0; it < refine_iters; ++it) {
    if (!isnan(civil_dawn_utc)) {
      double dec, EoT;
      solar_at_JD(JD0 + civil_dawn_utc/24.0, &dec, &EoT);
      double H = hour_angle_deg(lat_deg, dec, -6.0);
      double noon = 12.0 + (-lon_deg)*4.0/60.0 - EoT/60.0;
      civil_dawn_utc = norm24(noon - H/15.0);
    }
    if (!isnan(civil_dusk_utc)) {
      double dec, EoT;
      solar_at_JD(JD0 + civil_dusk_utc/24.0, &dec, &EoT);
      double H = hour_angle_deg(lat_deg, dec, -6.0);
      double noon = 12.0 + (-lon_deg)*4.0/60.0 - EoT/60.0;
      civil_dusk_utc = norm24(noon + H/15.0);
    }
    if (!isnan(sunrise_utc)) {
      double dec, EoT;
      solar_at_JD(JD0 + sunrise_utc/24.0, &dec, &EoT);
      double H = hour_angle_deg(lat_deg, dec, -0.833);
      double noon = 12.0 + (-lon_deg)*4.0/60.0 - EoT/60.0;
      sunrise_utc = norm24(noon - H/15.0);
    }
    if (!isnan(sunset_utc)) {
      double dec, EoT;
      solar_at_JD(JD0 + sunset_utc/24.0, &dec, &EoT);
      double H = hour_angle_deg(lat_deg, dec, -0.833);
      double noon = 12.0 + (-lon_deg)*4.0/60.0 - EoT/60.0;
      sunset_utc = norm24(noon + H/15.0);
    }
  }

  // 5) Populate outputs and flags
  out->have_civil_dawn = !isnan(civil_dawn_utc);
  out->have_civil_dusk = !isnan(civil_dusk_utc);
  out->have_sunrise    = !isnan(sunrise_utc);
  out->have_sunset     = !isnan(sunset_utc);

  out->polar_day   = (!out->have_sunrise && !out->have_sunset && out->have_civil_dawn);
  out->polar_night = (!out->have_sunrise && !out->have_sunset && !out->have_civil_dawn);

  out->civil_dawn_utc = civil_dawn_utc;
  out->civil_dusk_utc = civil_dusk_utc;
  out->sunrise_utc    = sunrise_utc;
  out->sunset_utc     = sunset_utc;

  // Local times
  {
    double tz = tz_offset_hours;
    out->civil_dawn_local = out->have_civil_dawn ? norm24(civil_dawn_utc + tz) : NAN;
    out->civil_dusk_local = out->have_civil_dusk ? norm24(civil_dusk_utc + tz) : NAN;
    out->sunrise_local    = out->have_sunrise    ? norm24(sunrise_utc + tz)    : NAN;
    out->sunset_local     = out->have_sunset     ? norm24(sunset_utc + tz)     : NAN;
  }

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <limits.h>


static void usage(const char *prog) {
  fprintf(stderr,
          "Usage: %s YYYY-MM-DD LAT LON TZ_OFFSET [LABEL]\n"
          "  LAT, LON: decimal degrees (N/E positive, S/W negative)\n"
          "  TZ_OFFSET: hours from UTC (e.g., -5, -4.0, 9)\n"
          "Example:\n"
          "  %s 2025-10-31 35.994 -78.898 -4 \"Durham NC\"\n",
          prog, prog);
}

static int parse_date(const char *s, int *y, int *m, int *d) {
  if (sscanf(s, "%d-%d-%d", y, m, d) != 3) return -1;
  if (*m < 1 || *m > 12) return -1;
  if (*d < 1 || *d > 31) return -1; /* fine to let compute_sun_times do strict checks */
  return 0;
}

static int parse_double(const char *s, double *out) {
  char *end = NULL;
  errno = 0;
  double v = strtod(s, &end);
  if (errno || end == s || *end != '\0') return -1;
  *out = v;
  return 0;
}

static int parse_int01(const char *s, int *out) {
  char *end = NULL;
  errno = 0;
  long v = strtol(s, &end, 10);
  if (errno || end == s || *end != '\0') return -1;
  if (v != 0 && v != 1) return -1;
  *out = (int)v;
  return 0;
}

// Pretty-printer used only in demo build
static void print_hhmm(double h){
  if (isnan(h)) { printf("--:--"); return; }
  int hh = (int)floor(h + 1e-12);
  int mm = (int)floor((h - hh)*60.0 + 0.5);
  if (mm == 60){ hh = (hh+1)%24; mm = 0; }
  printf("%02d:%02d \t", hh%24, mm);
}

// Demo CLI
int main(int argc, char **argv) {
  if (argc < 6 || (argc >= 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")))) {
    usage(argv[0]);
    return 2;
  }

  int year, month, day;
  double lat, lon, tz_off;
  int dst;

  if (parse_date(argv[1], &year, &month, &day) != 0) {
    fprintf(stderr, "Error: bad date '%s' (expected YYYY-MM-DD)\n", argv[1]);
    return 2;
  }
  if (parse_double(argv[2], &lat) != 0 || lat < -90.0 || lat > 90.0) {
    fprintf(stderr, "Error: bad latitude '%s' (range -90..90)\n", argv[2]);
    return 2;
  }
  if (parse_double(argv[3], &lon) != 0 || lon < -180.0 || lon > 180.0) {
    fprintf(stderr, "Error: bad longitude '%s' (range -180..180)\n", argv[3]);
    return 2;
  }
  if (parse_double(argv[4], &tz_off) != 0 || tz_off < -14.0 || tz_off > 14.0) {
    fprintf(stderr, "Error: bad TZ offset '%s' (hours from UTC, typically -12..+14)\n", argv[4]);
    return 2;
  }

  const char *label = (argc >= 6) ? argv[5] : NULL;

  SunTimes s;
  char err[160];

  int rc = compute_sun_times(year, month, day, lat, lon, tz_off, 1, &s, err, sizeof err);
  if (rc) {
    fprintf(stderr, "Error: %s\n", err);
    return rc;
  }

  printf("\nMico's Sun Times Calculator for the day of %04d-%02d-%02d\n", year, month, day);

  if (label) {
    printf("Location: %s\n", label);
  } else {
    printf("Location: (\xCF\x86=%.4f, \xCE\xBB=%.4f, TZ=%.2f, DST=%d)\n",
	   lat, lon, tz_off, dst);   // φ = 0xCF 0x86, λ = 0xCE 0xBB in UTF-8  
  }

  printf("\tDawn\tRise\tSet\tDusk\n");
  printf("Local \t");
  print_hhmm(s.civil_dawn_local); 
  print_hhmm(s.sunrise_local);    
  print_hhmm(s.sunset_local);     
  print_hhmm(s.civil_dusk_local);  printf("\n");

  printf("UTC\t");
  print_hhmm(s.civil_dawn_utc); 
  print_hhmm(s.sunrise_utc);    
  print_hhmm(s.sunset_utc);     
  print_hhmm(s.civil_dusk_utc);  printf("\n\n");
  
  if (s.polar_day)   printf("Note: Polar day (no true sunrise/sunset).\n");
  if (s.polar_night) printf("Note: Polar night.\n");

  return 0;
}
