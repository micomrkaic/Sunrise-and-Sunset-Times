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



function plotSunYearOctave()
% plotSunYearOctave — Octave/MATLAB compatible year plot with US/Eastern DST
% Requires: sunTimes.m (Octave-compatible version that returns serial datenums)
close all; clc;

% ===== USER SETTINGS =====
yearWanted = 2025;
lat_deg = 34.7128;
lon_deg = -84.0060;           % West negative
baseStdOffset = -5;           % US/Eastern standard offset hours
useUS_DST = true;             % apply US DST rules
% =========================

% Build LOCAL calendar midnights as serial datenumbers (no timezone objects)
localMidnights = datenum(yearWanted,1,1):1:datenum(yearWanted,12,31);
N = numel(localMidnights);

% Preallocate local-time event datenums (NaN means "no event")
dawnL    = NaN(N,1);
sunriseL = NaN(N,1);
sunsetL  = NaN(N,1);
duskL    = NaN(N,1);

polarDay   = false(N,1);
polarNight = false(N,1);

% DST boundaries (local time)
dstStart = us_dst_start_local(yearWanted);   % 02:00 on 2nd Sunday in March
dstEnd   = us_dst_end_local(yearWanted);     % 02:00 on 1st Sunday in November

% Iterate local calendar days
for k = 1:N
    LM = localMidnights(k);                                      % local midnight (datenum)
    % Offset at local midnight (hours)
    off_mid = local_offset_hours(LM, baseStdOffset, useUS_DST, dstStart, dstEnd);

    % Convert that local midnight instant to a UTC instant
    utcInstant = LM - off_mid/24;
    % Use the UTC *calendar date* (00:00 UTC) that contains this instant
    utcDate = floor(utcInstant);

    % Call sunTimes with tzOffsetHours=0 so outputs are UTC datenums
    out = sunTimesOctave(lat_deg, lon_deg, utcDate, 0);

    % Per-event conversion from UTC -> local using DST at the *event instant*
    if ~isnan(out.civilDawnUTC)
        off_evt = offset_at_utc_instant(out.civilDawnUTC, baseStdOffset, useUS_DST, dstStart, dstEnd);
        dawnL(k) = out.civilDawnUTC + off_evt/24;
    end
    if ~isnan(out.sunriseUTC)
        off_evt = offset_at_utc_instant(out.sunriseUTC, baseStdOffset, useUS_DST, dstStart, dstEnd);
        sunriseL(k) = out.sunriseUTC + off_evt/24;
    end
    if ~isnan(out.sunsetUTC)
        off_evt = offset_at_utc_instant(out.sunsetUTC, baseStdOffset, useUS_DST, dstStart, dstEnd);
        sunsetL(k) = out.sunsetUTC + off_evt/24;
    end
    if ~isnan(out.civilDuskUTC)
        off_evt = offset_at_utc_instant(out.civilDuskUTC, baseStdOffset, useUS_DST, dstStart, dstEnd);
        duskL(k) = out.civilDuskUTC + off_evt/24;
    end

    polarDay(k)   = out.flags.polarDay;
    polarNight(k) = out.flags.polarNight;

    % Ensure events belong to the same LOCAL day row (>= local midnight)
    if ~isnan(dawnL(k))    && dawnL(k)    < LM, dawnL(k)    = dawnL(k)    + 1.0; end
    if ~isnan(sunriseL(k)) && sunriseL(k) < LM, sunriseL(k) = sunriseL(k) + 1.0; end
    if ~isnan(sunsetL(k))  && sunsetL(k)  < LM, sunsetL(k)  = sunsetL(k)  + 1.0; end
    if ~isnan(duskL(k))    && duskL(k)    < LM, duskL(k)    = duskL(k)    + 1.0; end

    % Keep ordering within the same local solar day
    if ~isnan(dawnL(k))    && ~isnan(duskL(k))   && duskL(k)   < dawnL(k),   duskL(k)   = duskL(k)   + 1.0; end
    if ~isnan(sunriseL(k)) && ~isnan(sunsetL(k)) && sunsetL(k) < sunriseL(k), sunsetL(k) = sunsetL(k) + 1.0; end
end

% Convert events to time-of-day (hours since local midnight)
dawnH    = tod_hours(dawnL);
sunriseH = tod_hours(sunriseL);
sunsetH  = tod_hours(sunsetL);
duskH    = tod_hours(duskL);

% Stacked segments (hours)
night1 = dawnH;                                   % midnight -> dawn
civAM  = sunriseH - dawnH;                        % dawn -> sunrise
dayH   = sunsetH  - sunriseH;                     % sunrise -> sunset
civPM  = duskH    - sunsetH;                      % sunset -> dusk
night2 = 24 - duskH;                              % dusk -> midnight
segs = [night1, civAM, dayH, civPM, night2];
segs(segs < 0) = NaN;

% Polar visualization conventions
isAllNa = all(isnan(segs),2);
for k = 1:N
    if polarDay(k) && ~polarNight(k)
        segs(k,:) = [0 0 24 0 0];
    elseif polarNight(k) && ~polarDay(k)
        segs(k,:) = [24 0 0 0 0];
    elseif isAllNa(k)
        segs(k,:) = [NaN NaN NaN NaN NaN];
    end
end

%% === PLOT 1: Stacked areas (Night / Civil / Day / Civil / Night) ===
figure('Color','w','Name','Day Structure (Night/Twilight/Day)');
h = area(localMidnights, segs, 'LineStyle','none'); hold on;
if numel(h) == 5
    set(h(1),'FaceColor',[0.20 0.20 0.25]);   % Night
    set(h(2),'FaceColor',[1.00 0.70 0.30]);   % Civil AM
    set(h(3),'FaceColor',[1.00 0.90 0.35]);   % Day
    set(h(4),'FaceColor',[1.00 0.70 0.30]);   % Civil PM
    set(h(5),'FaceColor',[0.20 0.20 0.25]);   % Night
end
ylim([0 24]); set_hour_ticks(gca);
xlim([localMidnights(1) localMidnights(end)]);
grid on; box on;
title(sprintf('Sunlight Structure — %d  (lat %.4f, lon %.4f)', yearWanted, lat_deg, lon_deg));
ylabel('Local Time of Day');
legend({'Night','Civil Twilight (AM)','Day','Civil Twilight (PM)','Night'}, ...
       'Location','southoutside','Orientation','horizontal');
datetick('x','mmm','keeplimits');

%% === PLOT 2: Dots for Civil Dawn/Sunrise/Sunset/Civil Dusk ===
figure('Color','w','Name','Dawn/Sunrise/Sunset/Dusk (Dots)'); hold on;
p1 = plot(localMidnights, dawnH,    'LineStyle','none','Marker','o','MarkerSize',3);
p2 = plot(localMidnights, sunriseH, 'LineStyle','none','Marker','o','MarkerSize',3);
p3 = plot(localMidnights, sunsetH,  'LineStyle','none','Marker','o','MarkerSize',3);
p4 = plot(localMidnights, duskH,    'LineStyle','none','Marker','o','MarkerSize',3);
set(p1,'Color',[0.30 0.60 1.00]);   % dawn
set(p2,'Color',[0.00 0.45 0.74]);   % sunrise
set(p3,'Color',[0.85 0.33 0.10]);   % sunset
set(p4,'Color',[1.00 0.50 0.10]);   % dusk
ylim([0 24]); set_hour_ticks(gca);
xlim([localMidnights(1) localMidnights(end)]);
grid on; box on;
title(sprintf('Dawn / Sunrise / Sunset / Dusk — %d', yearWanted));
ylabel('Local Time of Day');
legend([p1 p2 p3 p4], {'Civil Dawn','Sunrise','Sunset','Civil Dusk'}, ...
       'Location','southoutside','Orientation','horizontal');
datetick('x','mmm','keeplimits');

drawnow;
end

%% ======================= LOCAL SUBFUNCTIONS =======================

function h = tod_hours(dnum)
% Time-of-day in hours from serial datenum (NaN-safe)
    h = 24 * (dnum - floor(dnum));
    h(isnan(dnum)) = NaN;
end

function set_hour_ticks(ax)
% Cross-version y-tick labels 00:00..24:00
    set(ax,'YLim',[0 24]);
    set(ax,'YTick',0:2:24);
    labs = arrayfun(@(x) sprintf('%02d:00', x), 0:2:24, 'uni', 0);
    set(ax,'YTickLabel',labs);
end

function off = local_offset_hours(local_dnum, baseStdOffset, useDST, dstStart, dstEnd)
% Offset (hours) at a given *local* time instant (datenum)
% US/Eastern: baseStdOffset = -5; DST adds +1 hour (i.e., -4)
    if ~useDST
        off = baseStdOffset; return;
    end
    if local_dnum >= dstStart && local_dnum < dstEnd
        off = baseStdOffset + 1;  % DST: -4
    else
        off = baseStdOffset;      % Standard: -5
    end
end

function off = offset_at_utc_instant(utc_dnum, baseStdOffset, useDST, dstStart, dstEnd)
% Determine local offset (hours) for a *UTC* instant.
% Iterate once: assume standard, convert to local, test DST, adjust.
    if ~useDST
        off = baseStdOffset; return;
    end
    local_guess = utc_dnum + baseStdOffset/24;
    if local_guess >= dstStart && local_guess < dstEnd
        off = baseStdOffset + 1;  % DST
    else
        off = baseStdOffset;      % Standard
    end
end

function dnum = us_dst_start_local(Y)
% Second Sunday in March, 02:00 local (US rules since 2007)
    d = datenum(Y,3,1);      % Mar 1
    dow = weekday(d);        % 1=Sun
    firstSun  = d + mod(1 - dow, 7);
    secondSun = firstSun + 7;
    dnum = secondSun + 2/24; % 02:00 local
end

function dnum = us_dst_end_local(Y)
% First Sunday in November, 02:00 local
    d = datenum(Y,11,1);     % Nov 1
    dow = weekday(d);
    firstSun = d + mod(1 - dow, 7);
    dnum = firstSun + 2/24;  % 02:00 local
end

