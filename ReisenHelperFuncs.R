CV = function(dat){ return(sd(dat)/mean(dat)) } #Calcs the coefficient of variation
YTrans = function(dat) { return(log(dat + 1)) } #Data transform selected by Reisen
YInvTrans = function(dat) { return(exp(dat)-1) } #Data backtransform selected by Reisen
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) } #Geometric mean
gm_mean2 = function(x, na.rm=TRUE){ exp(sum(log(x+1), na.rm=na.rm) / length(x))-1 } #Geometric mean
season = function(x) {
	x = as.numeric(x["WK"])
	if(x<23) {
		return("spring")
	} else if(x>35) {
		return("fall")
	} else {
		return("summer")
	}
}
CL = function(dat,alpha=0.05) {
	cent = mean(dat)
	ndat = length(dat)
	df = ndat-1
	lcl = cent + sd(dat)/sqrt(ndat)*qt(alpha/2,df)
	ucl = cent - sd(dat)/sqrt(ndat)*qt(alpha/2,df)
	return(paste(sprintf("%.2f",lcl),"-",sprintf("%.2f",ucl),sep=""))
}
tochr = function(x){return(sprintf("%02d",x))}
truncstr = function(x){return(substring(x,3))}
date_to_T = function(x) {
	d = as.Date(as.character(x), format="%Y%m%d")
	week = as.integer(format(d,"%U"))
	year = format(d,"%Y")
	ii = (year == "1995")
	week[ii] = week[ii] + 52
	samp = ceiling((week-16)/2) + 1
	return(samp)
}
yrwk_to_T = function(y,w) {
	d = as.Date(paste(y,w,0,sep=""), format="%Y%U%w")
	return(date_to_T(format(d,"%Y%m%d")))
}
dist_to_sea = function(lat,lon,bound) {
	lat1 = lat*pi/180
	lon1 = lon*pi/180
	lat2 = bound$Lat*pi/180
	lon2 = bound$Lon*pi/180
	delta_lat = lat2-lat1
	delta_lon = lon2-lon1
	a = sin(delta_lat/2)^2 + cos(lat1)*cos(lat2)*sin(delta_lon/2)^2
    c = 2*atan2(sqrt(a),sqrt(1-a))	
	dists = 6371 * c # in km
	return(min(dists))
}
