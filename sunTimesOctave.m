## /*
##  * This file is part of Mico's Sun Times Calculator
##  *
##  * Sun Times Calculator Is free software: you can redistribute it and/or modify
##  * it under the terms of the GNU General Public License as published by
##  * the Free Software Foundation, either version 3 of the License, or
##  * (at your option) any later version.
##  *
##  * Sun Times Calculator is distributed in the hope that it will be useful,
##  * but WITHOUT ANY WARRANTY; without even the implied warranty of
##  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##  * GNU General Public License for more details.
##  *
##  * You should have received a copy of the GNU General Public License
##  * along with Mico's Sun Times Calculator. If not, see <https://www.gnu.org/licenses/>.
##  */

function out = sunTimesOctave(lat_deg, lon_deg, date_ymd, tzOffsetHours)
% sunTimes  Compute civil dawn ("first light"), sunrise, civil dusk, and sunset.
%   out = sunTimes(lat_deg, lon_deg, date_ymd, tzOffsetHours)
%
% Compatibility: MATLAB R2013a+ and GNU Octave (no datetime/TimeZone needed)
%
% Inputs:
%   lat_deg        Latitude in degrees (+N).
%   lon_deg        Longitude in degrees (+E; West negative).
%   date_ymd       serial datenumber (datenum), or 'YYYY-MM-DD' string.
%                  (UTC calendar date used at 00:00)
%   tzOffsetHours  Local fixed offset from UTC in hours (e.g., -5). Default 0.
%
% Outputs (struct 'out'):  (all times are serial datenumbers; NaN means no event)
%   out.civilDawnUTC, out.sunriseUTC, out.sunsetUTC, out.civilDuskUTC
%   out.civilDawnLocal, out.sunriseLocal, out.sunsetLocal, out.civilDuskLocal
%   out.flags.polarDay, out.flags.polarNight   (logicals)
%
% Notes:
%   - Civil dawn/dusk use solar altitude h0 = -6 deg.
%   - Sunrise/sunset use h0 = -0.833 deg (refraction + Sun radius).
%   - Solar declination & EoT are evaluated at 00:00 UTC of the given UTC date
%     (typical error ≲ 1–2 minutes for rise/set).
%   - For DST, pass the appropriate offset for that date. (This function
%     intentionally uses a fixed numeric offset for Octave compatibility.)

    if nargin < 4, tzOffsetHours = 0; end

    % ---- Normalize input date to UTC midnight as serial datenumber ----
    if isnumeric(date_ymd) && isscalar(date_ymd)
        % Assume already a serial datenumber; force to start of that UTC day
        dateUTC = floor(date_ymd);
    elseif ischar(date_ymd) || (isstring(date_ymd) && isscalar(date_ymd))
        dateUTC = datenum(char(date_ymd), 'yyyy-mm-dd');  % 00:00 UTC that day
    else
        error('date_ymd must be datenum or ''YYYY-MM-DD'' string for Octave compatibility.');
    end

    % ---- Julian day / centuries since J2000.0 ----
    [Y,M,D] = datevec(dateUTC);                 % vectorized Y-M-D from datenum
    JD0 = gregorian2jd(Y,M,D);
    T0  = (JD0 - 2451545.0)/36525.0;

    % ---- Solar quantities (NOAA-style) ----
    M   = norm360(357.52911 + T0*(35999.05029 - 0.0001537*T0));              % mean anomaly (deg)
    C   = (1.914602 - T0*(0.004817 + 0.000014*T0))*sind(M) ...
        + (0.019993 - 0.000101*T0)*sind(2*M) + 0.000289*sind(3*M);           % equation of center (deg)
    L0  = norm360(280.46646 + T0*(36000.76983 + 0.0003032*T0));              % mean longitude (deg)
    lam = norm360(L0 + C);                                                    % ecliptic longitude (deg)
    eps = 23.439291 - T0*(0.013004167 + T0*(1.63889e-7 - 5.03611e-7*T0));    % obliquity (deg)

    % Declination
    dec = asind( sind(eps).*sind(lam) );

    % Equation of Time (minutes)
    eps_r = deg2rad(eps); L0_r = deg2rad(L0);
    e  = 0.016708634 - T0*(0.000042037 + 0.0000001267*T0);
    Mr = deg2rad(M);
    y  = tan(eps_r/2); y = y.*y;
    EoT_min = rad2deg( y.*sin(2*L0_r) - 2*e.*sin(Mr) + 4*e.*y.*sin(Mr).*cos(2*L0_r) ...
               - 0.5*y.*y.*sin(4*L0_r) - 1.25*e.*e.*sin(2*Mr) ) .* 4.0;

    % Local solar noon (UTC fractional hours)
    solarNoonUTC_hours = 12.0 + (-lon_deg)*4.0/60.0 - EoT_min/60.0;

    % Hour angles for dawn/dusk and sunrise/sunset
    H_civil = hourAngle(lat_deg, dec, -6.0);        % deg, NaN if none
    H_sun   = hourAngle(lat_deg, dec, -0.833);      % deg, NaN if none

    % ---- Build UTC times (as serial datenumbers; carry day offsets correctly) ----
    [civilDawnUTC, sunriseUTC, sunsetUTC, civilDuskUTC, polarDay, polarNight] = ...
        makeTimesUTC(dateUTC, solarNoonUTC_hours, H_civil, H_sun, lat_deg, dec);

    % ---- Shift to local via fixed offset (hours -> days) ----
    offDays = tzOffsetHours / 24.0;
    civilDawnLocal = civilDawnUTC + offDays;
    sunriseLocal   = sunriseUTC   + offDays;
    sunsetLocal    = sunsetUTC    + offDays;
    civilDuskLocal = civilDuskUTC + offDays;

    % ---- Ensure same local solar day ordering when both times exist ----
    if ~isnan(civilDawnLocal) && ~isnan(civilDuskLocal) && (civilDuskLocal < civilDawnLocal)
        civilDuskLocal = civilDuskLocal + 1.0;   % +1 day
    end
    if ~isnan(sunriseLocal) && ~isnan(sunsetLocal) && (sunsetLocal < sunriseLocal)
        sunsetLocal = sunsetLocal + 1.0;         % +1 day
    end

    % ---- Package outputs ----
    out = struct();
    out.civilDawnUTC   = civilDawnUTC;
    out.sunriseUTC     = sunriseUTC;
    out.sunsetUTC      = sunsetUTC;
    out.civilDuskUTC   = civilDuskUTC;

    out.civilDawnLocal = civilDawnLocal;
    out.sunriseLocal   = sunriseLocal;
    out.sunsetLocal    = sunsetLocal;
    out.civilDuskLocal = civilDuskLocal;

    out.flags = struct('polarDay',logical(polarDay),'polarNight',logical(polarNight));
end

% ---------- helpers ----------

function JD = gregorian2jd(Y,M,D)
    a = floor((14 - M)/12);
    y = Y + 4800 - a;
    m = M + 12*a - 3;
    JD = D + floor((153*m + 2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045;
end

function x = norm360(x)
    x = mod(x,360); if x < 0, x = x + 360; end
end

function Hdeg = hourAngle(lat_deg, dec_deg, alt_deg)
% Returns hour angle (deg) in [0,180] at which solar altitude equals alt_deg.
% Returns NaN if never reaches that altitude on this date (polar conditions).
    lat = deg2rad(lat_deg);
    dec = deg2rad(dec_deg);
    alt = deg2rad(alt_deg);
    cosH = (sin(alt) - sin(lat).*sin(dec)) ./ (cos(lat).*cos(dec));

    % Decide "no solution" vs small numeric overshoot; then clamp
    eps = 1e-12;
    if any(isnan(cosH)) || cosH > 1 + eps || cosH < -1 - eps
        Hdeg = NaN; return
    end
    cosH = min(1,max(-1,cosH));
    Hdeg = acosd(cosH);   % [0..180]
end

function [civilDawnUTC, sunriseUTC, sunsetUTC, civilDuskUTC, polarDay, polarNight] = ...
         makeTimesUTC(dateUTC_dnum, solarNoonHr, Hcivil, Hsun, lat_deg, dec_deg)

    civilDawnUTC = NaN; sunriseUTC = NaN; sunsetUTC = NaN; civilDuskUTC = NaN;
    polarDay = false; polarNight = false;

    % Civil twilight times (if exist)
    if ~isnan(Hcivil)
        civilDawnUTC = fracHourToDatenum(dateUTC_dnum, solarNoonHr - Hcivil/15);
        civilDuskUTC = fracHourToDatenum(dateUTC_dnum, solarNoonHr + Hcivil/15);
    end

    % Sunrise/Sunset times (if exist)
    if ~isnan(Hsun)
        sunriseUTC = fracHourToDatenum(dateUTC_dnum, solarNoonHr - Hsun/15);
        sunsetUTC  = fracHourToDatenum(dateUTC_dnum, solarNoonHr + Hsun/15);
    else
        % No rise/set today: classify
        if isnan(Hcivil)
            % Neither civil twilight nor rise/set -> deep polar condition.
            % Use noon altitude to decide day vs night.
            hNoon = noonAltitude(lat_deg, dec_deg);  % degrees
            if hNoon > -0.833
                polarDay = true;   % always above -0.833 (effectively above horizon all day)
            else
                polarNight = true; % always below -0.833 and even below -6 at noon -> deep night
            end
        else
            % Civil twilight exists but no rise/set -> polar night with twilight.
            polarNight = true;
        end
    end
end

function dnum = fracHourToDatenum(dateUTC_dnum, hourFrac)
% Convert fractional UTC hours from the *given UTC date's* midnight (serial datenum),
% carrying day offsets correctly (works for hourFrac < 0 or >= 24).
    dOff  = floor(hourFrac/24);        % whole-day offset (can be negative)
    hRest = hourFrac - 24*dOff;        % remainder within [0,24)
    dnum  = floor(dateUTC_dnum) + dOff + hRest/24;
end

function hdeg = noonAltitude(lat_deg, dec_deg)
% Solar altitude (deg) at local solar noon (hour angle = 0)
    lat = deg2rad(lat_deg);
    dec = deg2rad(dec_deg);
    h = asin( sin(lat).*sin(dec) + cos(lat).*cos(dec) ); % radians
    hdeg = rad2deg(h);
end

