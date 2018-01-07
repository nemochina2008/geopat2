/****************************************************************************
 *
 * MODULE:	p.sim.segment2
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski

 * PURPOSE:	information retrival using categorical maps:
 *		compares grid of histograms
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#define MAIN
#include "local_proto.h"
#include <grass/glocale.h>

int main(int argc, char *argv[])
{
	struct GModule *module;

	/* menu */
	struct Option* opt_input_grids=add_menu_item(INPUT_GRIDS);
	struct Option* opt_nulls=add_menu_item(OPTION_NULLS);
	struct Option* opt_dist_matrix_file=add_menu_item(INPUT_DIST_MATRIX_FILE);
	struct Option* opt_measure=add_menu_item(OPTION_MEASURE);
	struct Option* opt_lthreshold=add_local_menu_item(OPTION_LTHRESHOLD);
	struct Option* opt_uthreshold=add_local_menu_item(OPTION_UTHRESHOLD);
	struct Option* opt_weights=add_local_menu_item(OPTION_WEIGHTS);
	struct Option* opt_swap=add_local_menu_item(OPTION_SWAP);
	struct Option* opt_refine=add_local_menu_item(OPTION_REFINE);
	struct Option* opt_minarea=add_menu_item(OPTION_MINAREA);
	struct Option* opt_maxhist=add_menu_item(OPTION_MAXHIST);
	struct Option* opt_output_layer=add_menu_item(OUTPUT_LAYER);

	struct Flag* flag_complete=add_local_menu_flag(FLAG_COMPLETE);
	struct Flag* flag_threshold=add_local_menu_flag(FLAG_THRESHOLD);
	struct Flag* flag_growing=add_local_menu_flag(FLAG_GROWING);
	struct Flag* flag_hierarhical=add_local_menu_flag(FLAG_HIERARHICAL);
	struct Flag* flag_all=add_local_menu_flag(FLAG_ALL);
	struct Flag* flag_quad=add_local_menu_flag(FLAG_QUAD);

	struct Cell_head window;
	CELL* results=NULL;
	struct area** areas=NULL;
	HEXGRID* hexgrid;
	DATAINFO** datainfo;
	LOCAL_PARAMS* parameters;
	parameters=G_malloc(sizeof(LOCAL_PARAMS));
	datainfo=G_malloc(sizeof(DATAINFO));

	G_gisinit(argv[0]);
	module = G_define_module();
	G_add_keyword(_("similarity"));
	G_add_keyword(_("segmentation"));
	G_add_keyword(_("information retrieval"));
    module->description =
	_("segments grid of histograms using available similarity measures and create layer of unique regions");
	if (G_parser(argc, argv))
		exit(EXIT_FAILURE);

	if(flag_growing->answer && flag_hierarhical->answer)
		G_fatal_error(_("Only one flag -h or -g can be used"));

	if(strcmp(opt_measure->answer,"emd")==0 && opt_dist_matrix_file->answer==NULL)
		G_fatal_error(_("Eart Mover Distance (emd) needs a distances matrix. Add the parameter: distfile!"));

	parameters->null_threshold=atof(opt_nulls->answer);
	if(parameters->null_threshold<0 || parameters->null_threshold>1)
		G_fatal_error(_("Null threshold must be between 0 and 1"));
	if(parameters->null_threshold>0.75)
		G_warning("Null threshold <%f> may be too big",parameters->null_threshold);

	parameters->lower_similarity_threshold=atof(opt_lthreshold->answer);
	if(parameters->lower_similarity_threshold<0 || parameters->lower_similarity_threshold>1)
		G_fatal_error(_("Lower distance threshold must be between 0 and 1"));

	parameters->upper_similarity_threshold=atof(opt_uthreshold->answer);
	if(parameters->upper_similarity_threshold<0 || parameters->upper_similarity_threshold>1)
		G_fatal_error(_("Upper distance threshold must be between 0 and 1"));

	if(parameters->upper_similarity_threshold < parameters->lower_similarity_threshold)
		G_fatal_error(_("Upper distance threshold cannot be smaller than lower threshold"));

	parameters->sampling_threshold=atoi(opt_maxhist->answer);
	if(parameters->sampling_threshold<0)
		G_fatal_error(_("Sampling threshold cannot be negative"));

	parameters->minarea=atoi(opt_minarea->answer);
	if(parameters->minarea<0)
		G_fatal_error(_("area must be non-negative"));

	parameters->refinement=atoi(opt_refine->answer);

	parameters->reduction=0;

	parameters->swap_threshold=atof(opt_swap->answer);
	if(parameters->swap_threshold>1)
		G_fatal_error(_("swap must be between 0 and 1 or negative to skip"));

	parameters->quad_mode=(flag_quad->answer!=0);
	parameters->complete_linkage=(flag_complete->answer!=0);
	parameters->all_layers=(flag_all->answer!=0);
	if(parameters->all_layers && opt_weights->answer)
		G_warning("ignore weigths in <all layers> mode");

	if(parameters->quad_mode && parameters->refinement>=0) {
		G_warning("ignore refinement in quad model");
		parameters->refinement=0;
	}


	Rast_get_window(&window);
	parameters->calculate=measure_function(opt_measure->answer);

	int i;
	int num_of_layers=get_num_of_grids(opt_input_grids->answers);
	datainfo=malloc(num_of_layers*sizeof(DATAINFO*));
	for(i=0;i<num_of_layers;++i) {
		datainfo[i]=malloc(num_of_layers*sizeof(DATAINFO));
		init_grid_datainfo(datainfo[i],opt_input_grids->answers[i],"OUTPUT",0);
		if(i)
			compare_grids_datainfo(datainfo[i-1], datainfo[i]);
	}


	if(strcmp(opt_measure->answer,"emd")==0) {
		if(!measure_parameters_dist_matrix(parameters->parameters,opt_dist_matrix_file->answer))
			G_fatal_error(_("Problem with distance matrix reading"));
	} else if(strcmp(opt_measure->answer,"dtw")==0) {
		measure_parameters_vector_size(parameters->parameters,datainfo[0]->pattern_size);
	}

	G_message("Read data...");

	for(i=0;i<num_of_layers;++i)
		read_histograms_to_memory(datainfo[i],parameters);

	hexgrid=hex_build_topology(datainfo,parameters,num_of_layers,opt_weights->answer,parameters->quad_mode);
	areas=hex_build_areas(datainfo,hexgrid,parameters);
	results=hex_init_results(hexgrid);
	parameters->parameters=init_measure_parameters(datainfo[0]->size_of_histogram,0); /* we will use distance instead of similarity */

	/* seeding starts here */
	unsigned num_of_seeds;
	unsigned* seeds=hex_find_seeds(hexgrid,parameters,areas,&num_of_seeds);

	/* thresholds only */
	if(flag_threshold->answer!=0) {
		double* thresholds=create_thresholds_map(datainfo[0],hexgrid,parameters,areas);
		Rast_get_window(&window);
		G_message(_("Change window to write results"));
		Rast_set_window(&(datainfo[0]->cell_hd));
		char map_name[100]="\0";
		sprintf(map_name,"%s_threshold",opt_output_layer->answer);
		write_map(datainfo[0],parameters,(void*)thresholds,map_name,DCELL_TYPE,2);
		G_message(_("Restore original region definition"));
		Rast_set_window(&window);
		free(seeds);
		hex_remove_hexgrid(hexgrid);
		free(results);
		free(datainfo);
		free(parameters); /* TODO: ADD RELEASE FUNCTION */;
		free(hexgrid);
		exit(EXIT_SUCCESS);
	}

	/* segmentation starts here */
	if(!flag_growing->answer)
		hex_region_growing(hexgrid,parameters,areas,results,seeds,num_of_seeds);

	if(!flag_hierarhical->answer)
		hex_hierarhical(hexgrid,parameters,areas,results);

	if(parameters->minarea>0)
		hex_minarea(hexgrid,parameters,areas,results);

	if(parameters->swap_threshold>=0) {
		swap_areas(hexgrid,parameters,areas,results,0);
		if(parameters->minarea>0)
			hex_minarea(hexgrid,parameters,areas,results);
	}

	/* refinement */
	if (parameters->refinement>=0) {
		G_message(_("Refinement..."));
		free(results);
		results=hex_create_segment_map(datainfo[0],hexgrid,parameters,areas); /* rebuild result in quad topology */
		hex_remove_hexgrid(hexgrid); /* clean hexgrid */
		/* rebuild input in quad topology */
		hexgrid=hex_build_topology(datainfo,parameters,num_of_layers,opt_weights->answer,1);
		areas=hex_rebuild_areas(hexgrid,parameters,results);
		swap_areas(hexgrid,parameters,areas,results,parameters->refinement);
		parameters->minarea*=2; /* temporary fixed */
		parameters->upper_similarity_threshold=1;
		if(parameters->minarea>0)
			hex_minarea(hexgrid,parameters,areas,results);
	}

	hex_reclass(hexgrid,areas);
	int* segment_map=hex_create_segment_map(datainfo[0],hexgrid,parameters,areas);


	Rast_get_window(&window);
	G_message(_("Change window to write results"));
	Rast_set_window(&(datainfo[0]->cell_hd));
	char map_name[100]="\0";
	sprintf(map_name,"%s_SEGMENT",opt_output_layer->answer);
	write_map(datainfo[0],parameters,(void*)segment_map,map_name,CELL_TYPE,1); /* write segment map */
	G_message(_("Restore original region definition"));
	Rast_set_window(&window);

/*	for(i=0;i<ncells;++i)
		if(areas[i]!=NULL)
			remove_area(areas+i);
	free(areas);
*/
	hex_remove_hexgrid(hexgrid);
	free(segment_map);
	free(results);
	free(datainfo);
	free(parameters); /* TODO: ADD RELEASE FUNCTION */;
	free(hexgrid);


exit(EXIT_SUCCESS);
}

