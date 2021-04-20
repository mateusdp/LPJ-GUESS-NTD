//
//  fluxnet.h
//  guess
//
//  Created by Adrian Gustafson on 2018-06-14.
//
//

#ifndef LPJ_GUESS_FLUXNET_H
#define LPJ_GUESS_FLUXNET_H

#include "guess.h"
#include "inputmodule.h"
#include <vector>
#include "gutil.h"
#include "globalco2file.h"
#include "soilinput.h"
#include "spinupdata.h"
#include "cru_ts30.h"
#include "lamarquendep.h"
#include "externalinput.h"

class FluxnetInput : public InputModule {
public:
	FluxnetInput();

	~FluxnetInput();

    void init();
    
    bool getgridcell(Gridcell& gridcell);
    
    bool getclimate(Gridcell& gridcell);


	/// See base class for documentation about this function's responsibilities
	void getlandcover(Gridcell& gridcell);

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) { management_input.getmanagement(gridcell); }

	// Constants associated with historical climate data set

	/// number of years of historical climate
	static const int NYEAR_HIST = CRU_FastArchive::NYEAR_HIST;

	/// calendar year corresponding to first year in data set
	static const int FIRSTHISTYEAR = CRU_FastArchive::FIRSTHISTYEAR;

	/// number of years to use for temperature-detrended spinup data set
	/** (not to be confused with the number of years to spinup model for, which
	* is read from the ins file)
	*/
	static const int NYEAR_SPINUP_DATA = 30;
	/*
	void adjust_raw_forcing_data(double hist_mtemp[NYEAR_HIST][12], 
								 double hist_mprec[NYEAR_HIST][12],
								 double hist_msun[NYEAR_HIST][12], 
								 double fluxnet_temp[12],
								 double fluxnet_prec[12],
								 double fluxnet_sun[12]);
    */
protected:

	/// Gets monthly ndep values for a given calendar year
	/** To be used by sub-classes that wish to do their own
	*  distribution of monthly values to daily values.
	*
	*  The ndep values returned are for the current gridcell,
	*  i.e. the gridcell chosen in the most recent call to
	*  getgridcell().
	*
	*  \param calendar_year The calendar (not simulation!) year for which to get ndep
	*  \param mndrydep      Pointer to array holding 12 doubles
	*  \param mnwetdep      Pointer to array holding 12 doubles
	*/
	void get_monthly_ndep(int calendar_year,
	                      double* mNHxdrydep, double* mNOydrydep,
						  double* mNHxwetdep, double* mNOywetdep);

	/// Gives sub-classes a chance to modify the forcing data
	/** This function will be called just after the forcing data for the historical
	*  period has been read in for a gridcell. Sub-classes can override this function
	*  and modify the data if needed, for instance adjusting according to site data.
	*
	*  Note that modifying this data will also affect the spinup period since
	*  the spinup forcing is based on the historical period.
	*
	*  \param hist_mtemp  Monthly temperature values for each year
	*  \param hist_mprec  Monthly precipitation values for each year
	*  \param hist_msun   Monthly sunshine values for each year
	*/
	void adjust_raw_forcing_data(double hist_mtemp[NYEAR_HIST][12],
		double hist_mprec[NYEAR_HIST][12],
		double hist_msun[NYEAR_HIST][12],
		double fluxnet_temp[12],
		double fluxnet_prec[12],
		double fluxnet_sun[12]);


private:
	/// Type for storing grid cell longitude, latitude and description text
	struct Coord {

		int id;
		double lon;
		double lat;
		xtring descrip;
	};

	std::vector<std::pair<double, double> > translate_gridlist_to_coord(ListArray_id<Coord>& gridlist);

	// Soil input module
	SoilInput soilinput;
	/// Land cover input module
	LandcoverInput landcover_input;
	/// Management input module
	ManagementInput management_input;

	/// search radius to use when finding CRU data
	double searchradius;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// Flag for getgridcell(). True indicates that the first gridcell has not been read yet by getgridcell()
	bool first_call;

	// Timers for keeping track of progress through the simulation
	Timer tprogress, tmute;
	static const int MUTESEC = 20; // minimum number of sec to wait between progress messages

								   /// Yearly CO2 data read from file
								   /**
								   * This object is indexed with calendar years, so to get co2 value for
								   * year 1990, use co2[1990]. See documentation for GlobalCO2File for
								   * more information.
								   */
	GlobalCO2File co2;

	/// Monthly temperature for current grid cell and historical period
	double hist_mtemp[NYEAR_HIST][12];

	/// Monthly precipitation for current grid cell and historical period
	double hist_mprec[NYEAR_HIST][12];

	/// Monthly sunshine for current grid cell and historical period
	double hist_msun[NYEAR_HIST][12];

	/// Monthly frost days for current grid cell and historical period
	double hist_mfrs[NYEAR_HIST][12];

	/// Monthly precipitation days for current grid cell and historical period
	double hist_mwet[NYEAR_HIST][12];

	/// Monthly DTR (diurnal temperature range) for current grid cell and historical period
	double hist_mdtr[NYEAR_HIST][12];

	/// Monthly mean wind for current grid cell and historical period
	double hist_mwind[NYEAR_HIST][12];

	/// Monthly mean relative humidity for current grid cell and historical period
	double hist_mrhum[NYEAR_HIST][12];

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Spinup data for current grid cell - temperature
	Spinup_data spinup_mtemp;
	/// Spinup data for current grid cell - precipitation
	Spinup_data spinup_mprec;
	/// Spinup data for current grid cell - sunshine
	Spinup_data spinup_msun;

	/// Spinup data for current grid cell - frost days
	Spinup_data spinup_mfrs;
	/// Spinup data for current grid cell - precipitation days
	Spinup_data spinup_mwet;
	/// Spinup data for current grid cell - DTR (diurnal temperature range)
	Spinup_data spinup_mdtr;
	/// Spinup data for current grid cell - wind
	Spinup_data spinup_mwind;
	/// Spinup data for current grid cell - relative humidity
	Spinup_data spinup_mrhum;

	/// Daily temperature for current year
	double dtemp[Date::MAX_YEAR_LENGTH];
	/// Daily precipitation for current year
	double dprec[Date::MAX_YEAR_LENGTH];
	/// Daily sunshine for current year
	double dsun[Date::MAX_YEAR_LENGTH];
	// Daily diurnal temperature range for current year
	double ddtr[Date::MAX_YEAR_LENGTH];
	/// Daily N deposition for current year
	double dndep[Date::MAX_YEAR_LENGTH];
	// Daile mean wind
	double dwind[Date::MAX_YEAR_LENGTH];
	// Daily mean relative humidity
	double drhum[Date::MAX_YEAR_LENGTH];
	/// Daily N deposition for current year
	double dNH4dep[Date::MAX_YEAR_LENGTH];
	double dNO3dep[Date::MAX_YEAR_LENGTH];

};

#endif // LPJ_GUESS_FLUXNET_H
