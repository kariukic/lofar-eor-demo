#!/usr/bin/env nextflow

include {
    getTime
} from './utils.nf'

process setupDirs {
    debug true

    input:
        val ready

    output:
        val true

    script:

        """
            mkdir -p ${params.data.path}/${params.out.logs}/flag
            mkdir -p ${params.data.path}/${params.out.logs}/aoquality
            mkdir -p ${params.data.path}/${params.out.logs}/calibrate
            mkdir -p ${params.data.path}/${params.out.logs}/average
            mkdir -p ${params.data.path}/${params.out.logs}/image
            mkdir -p ${params.data.path}/${params.out.logs}/datalists
            mkdir -p "${params.data.path}/${params.out.logs}/power_spectrum/di"
            mkdir -p "${params.data.path}/${params.out.logs}/power_spectrum/bp"
            mkdir -p "${params.data.path}/${params.out.logs}/power_spectrum/dd"
            mkdir -p "${params.data.path}/${params.out.logs}/power_spectrum/pdd"
            mkdir -p ${params.data.path}/${params.out.results}
        """
}


process FlagIntra {
    debug true
    label 'sing'

    input:
        val ready
        path ms

    output:
        val true

    script:
        time = getTime()
        """
        python3 /home/codex/chege/software/pipelines/twoleap/templates/flag_intrastations.py -i ${ms} > "${params.data.path}/${params.out.logs}/flag/flag_intrastations_${ms}_${time}.log" 2>&1
        """

}

process FilterInter {
    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
        val ready
        path ms

    output:
        path "${ms.getName() + '.noInter'}"

    script:
        time = getTime()
        """
        DP3 msin=${ms} steps=[filter] filter.remove=true filter.baseline="[CR]S*&&" msout="${ms.getName() + '.noInter'}" msout.overwrite=True > "${params.data.path}/${params.out.logs}/flag/filter_intrastations_and_international_stations_${ms}_${time}.log" 2>&1
        """
}


process AOFlag {

    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
        val ready
        path ms
        val column
        path aoflagger_strategy
        val interpolate

    output:
        val true , emit: done
    

    script:
        time=getTime()

        if ( interpolate == 1 )
            """
            DP3 steps=[aoflag,interpolate] msin=${ms} msin.datacolumn=${column} aoflag.type=aoflagger aoflag.strategy=${aoflagger_strategy} msout=. msout.overwrite=True > "${params.data.path}/${params.out.logs}/flag/flag_ao_${ms}_${column}_${time}.log" 2>&1
            """
        else
            """
            DP3 steps=[aoflag] msin=${ms} msin.datacolumn=${column} aoflag.type=aoflagger aoflag.strategy=${aoflagger_strategy} msout=. msout.overwrite=True > "${params.data.path}/${params.out.logs}/flag/flag_ao_${ms}_${column}_${time}.log" 2>&1
            """
}

process AverageData {
    label 'sing'
    publishDir params.data.path, mode: 'move'

    input:
        val ready
        tuple path(msin), val(msout)
        val data_column
        val timestep
        val freqstep


    output:
        path "${msout}", emit: averaged_ms
        val true , emit: done_averaging

    script:
        time=getTime()
        """
        DP3 steps=[avg] msin=${msin} msin.datacolumn=${data_column} msout=${msout} avg.type=average avg.timestep=${timestep} avg.freqstep=${freqstep} msout.overwrite=True > "${params.data.path}/${params.out.logs}/average/average_${msin}_${timestep}tstep_${freqstep}freqstep_${time}.log" 2>&1
        """
}


process WScleanImage {
    label 'sing'
    publishDir "${params.data.path}/${params.out.results}/wsclean/${datacol}", pattern: "*.fits", mode: "move", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/wsclean/${datacol}/plots", pattern: "*.png", mode: "move", overwrite: true

    // publishDir "${params.data.path}/${params.out.results}/images", pattern: "*.txt", mode: "copy", overwrite: true

    input:
        val ready
        tuple val(mses), val(imname)
        val size
        val scale
        val niter
        val pol
        val chansout
        val minuvl
        val maxuvl
        val weight
        val polfit
        val datacol
        

    output:
        path "*.fits"
        path "*.png"
        val true , emit: done
        // path "${imname}-sources.txt", emit: model

    script:
        if ( chansout == 1 )
            """
            wsclean -v -log-time -name ${imname} -data-column ${datacol} -pol ${pol} -weight ${weight} -scale ${scale} -size ${size} ${size} -make-psf -niter ${niter} -gridder wgridder -reorder ${mses} > ${params.data.path}/${params.out.logs}/image/wsclean_${imname}_image.log 2>&1

            python3 ${projectDir}/templates/plot_images.py plot --images '*-I-image.fits' --filename ${imname}
            """

        else
            """
            wsclean -v -log-time -name ${imname} -data-column ${datacol} -pol ${pol} -weight ${weight} -scale ${scale} -size ${size} ${size} -niter ${niter} -apply-primary-beam -make-psf -join-channels -channels-out ${chansout} -gridder wgridder -no-update-model-required -no-dirty -no-mf-weighting ${mses} > ${params.data.path}/${params.out.logs}/image/wsclean_${imname}_image.log 2>&1

            ls *-I-image.fits > imlist.txt

            python3 ${projectDir}/templates/plot_images.py plot --imagelist imlist.txt --filename ${imname}
            """

}


process AOqualityCollect {
    label 'sing'

    input:
        val ready
        path full_ms_path
        val data_column

    output:
        val true

    script:
        time=getTime()
        """
        aoquality collect -d ${data_column} ${full_ms_path}  > ${params.data.path}/${params.out.logs}/aoquality/collect_${full_ms_path.getName()}_${time}.log 2>&1
        """
}

process AOqualityCombine {
    label 'sing'
    publishDir "${params.data.path}/${params.out.results}/aoquality/${output_name}", pattern: "*.qs", mode: "copy", overwrite: true
    publishDir "${params.data.path}/${params.out.results}/aoquality/${output_name}/plots", pattern: "*.pdf", mode: "copy", overwrite: true
    
    input:
        val ready
        val file_list
        val output_name

    output:
        path "${output_name}.qs"
        path "*.pdf"
        val true, emit: done

    script:
        time=getTime()
        List tlist = file(file_list).readLines()
        String mses = tlist.collect {"${it}"}.join(" ")
        """
        aoquality combine ${output_name}.qs ${mses} > ${params.data.path}/${params.out.logs}/aoquality/combine_${output_name}_${time}.log 2>&1
        python3 ${projectDir}/templates/plot_flags.py plot_occ ${file_list} --filename ${output_name} >> ${params.data.path}/${params.out.logs}/aoquality/combine_${output_name}_${time}.log 2>&1
        python3 ${projectDir}/templates/plot_aoqstats.py plot_aoq "${output_name}.qs" --name ${output_name} >> ${params.data.path}/${params.out.logs}/aoquality/combine_${output_name}_${time}.log 2>&1
        """
}


process GetData {
    debug true

    input:
        val ready
        val glob
        val txtname

    output:
        path "${txtname}"

    script:
        """
        ls -d ${glob} >> ${txtname}
        cp ${txtname} ${params.data.path}
        """
}

// process ReadMSlist {

//     input:
//         val ready
//         val txt

//     output:
//         val mslist
    
//     script:
//         mslist = readTxtIntoString("${txt}")
//         """
//         """
// }

// def readTxtIntoString (txt) {
//     List tlist = file(txt).readLines()
//     String tstring = tlist.collect {"${it}"}.join(" ")

//     return tstring
// }


process DICalibrate {
    debug true
    label 'sing'
    maxForks 5
    publishDir "${ms}" , mode: 'copy'

    input:
        val ready
        path ms
        path parset
        path sourcedb
        val solsfile
        val incol
        val solint
        val uvlambdamin
        val uvlambdamax
        val nchan

    output:
        path "${solsfile}"

    script:

        time = getTime()

        """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.nchan=${nchan} > "${params.data.path}/${params.out.logs}/calibrate/di_cal_${ms}_${solsfile}_${time}.log" 2>&1
        """
}


process ApplyGains {
    debug true
    label 'sing'

    input:
        val ready
        tuple path(ms), path(solsfile)
        val parset
        val incol
        val outcol

    output:
        path "${ms}"

    script:
        
        time = getTime()

        """
        DP3 ${parset} msin=${ms} applycal.parmdb=${solsfile} msin.datacolumn=${incol} msout.datacolumn=${outcol} > "${params.data.path}/${params.out.logs}/calibrate/applygains_di_${ms}_${solsfile}_to_${outcol}_${time}.log" 2>&1
        """
}

process H5ParmCollect {
    debug true
    label 'sing'

    publishDir "${params.data.path}/results/solutions/${output_name}", pattern: "*.h5", mode: "move", overwrite: true
    publishDir "${params.data.path}/results/solutions/${output_name}/plots", pattern: "*.png", mode: "move", overwrite: true

    input:
        val ready
        val solution_files
        val output_name

    output:
        path "${output_name}.h5", emit: combined_sols
        path "*.png", emit: plots
        val true, emit: done

    script:
        """
        H5parm_collector.py ${solution_files} -o ${output_name}.h5 > "${params.data.path}/${params.out.logs}/calibrate/h5parm_collect.log" 2>&1
        #soltool plot --plot_dir \$(pwd) ${output_name}.h5 >> "${params.data.path}/${params.out.logs}/calibrate/h5parm_collect.log" 2>&1
        python3 ${projectDir}/templates/plot_sols.py plot ${output_name}.h5 >> "${params.data.path}/${params.out.logs}/calibrate/h5parm_collect.log" 2>&1
        """
}


process DIGainCal {
    debug true
    label 'sing'
    maxForks 5
    publishDir "${ms}" , mode: 'copy'

    input:
        val ready
        path ms
        path parset
        path sourcedb
        val solsfile
        val incol
        val solint
        val uvlambdamin
        val uvlambdamax
        val nchan

    output:
        path "${solsfile}"

    script:

        time = getTime()

        """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} gaincal.sourcedb=${sourcedb} gaincal.parmdb=${solsfile} gaincal.solint=${solint} gaincal.uvlambdamin=${uvlambdamin} gaincal.uvlambdamax=${uvlambdamax} gaincal.nchan=${nchan} > "${params.data.path}/${params.out.logs}/calibrate/gain_cal_di_${ms}_${solsfile}_${time}.log"
        """
}


process DDCalibrate {
    debug true
    label 'sing'
    maxForks 5
    publishDir "${ms}" , mode: 'copy'

    input:
        val ready
        path ms
        path parset
        path sourcedb
        val solsfile
        val incol
        val calmode
        val solint
        val uvlambdamin
        val uvlambdamax
        val nchan
        val usebeam
        val beammode
        val smoothnessconstraint
        val truncateksmoothkernel
        val robust_reg
        val propagate_sols
        val maxiter
        val beamproximitylimit
        val correctfreqsmearing
        val flagstations

    output:
        path "${solsfile}"

    script:

        time = getTime()

        if ( flagstations )

            """
            DP3 ${parset} steps=[preflagger,ddecal] msin=${ms} preflagger.baseline="${flagstations}" msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.model_weighted_constraints=${robust_reg} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.correctfreqsmearing=${correctfreqsmearing} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} > "${params.data.path}/${params.out.logs}/calibrate/cal_dd_${ms}_${solsfile}_${time}.log"
            """

        else

            """
            DP3 ${parset} steps=[ddecal] msin=${ms} msin.datacolumn=${incol} ddecal.sourcedb=${sourcedb} ddecal.h5parm=${solsfile} ddecal.solint=${solint} ddecal.uvlambdamin=${uvlambdamin} ddecal.uvlambdamax=${uvlambdamax} ddecal.nchan=${nchan} ddecal.mode=${calmode} ddecal.smoothnessconstraint=${smoothnessconstraint} ddecal.model_weighted_constraints=${robust_reg} ddecal.propagatesolutions=${propagate_sols} ddecal.maxiter=${maxiter} ddecal.beamproximitylimit=${beamproximitylimit} ddecal.correctfreqsmearing=${correctfreqsmearing} ddecal.usebeammodel=${usebeam} ddecal.beammode=${beammode} ddecal.smoothness_kernel_truncation=${truncateksmoothkernel} > "${params.data.path}/${params.out.logs}/calibrate/cal_dd_${ms}_${solsfile}_${time}.log"
            """
}


process SubtractSources {
    label 'sing'
    input:
        val ready
        tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file) //, path(sources_to_subtract_file)
        path subtraction_parset
        val input_datacolumn
        val output_datacolumn

    output:
        path "${full_ms_path}"
    
    script:
        time = getTime()

        """
        DP3 ${subtraction_parset} msin=${full_ms_path} sub.applycal.parmdb=${calibration_solutions_file} sub.sourcedb=${sourcedb_name} msin.datacolumn=${input_datacolumn} msout.datacolumn=${output_datacolumn} > "${params.data.path}/${params.out.logs}/calibrate/subtract_dd_${full_ms_path}_dd_subtract_${time}.log" 2>&1
        """
}


process UVWFlag {

    debug true
    label 'sing'
    publishDir "${params.data.path}", mode: 'move'

    input:
        val ready
        path msin
        val data_column
        val uvlambdamin
        val uvlambdamax

    output:
        path "${msin}.l${uvlambdamin}to${uvlambdamax}", emit: msout
        val true , emit: done

    script:
        time=getTime()
        """
        DP3 steps=[uvwflag] msin=${msin} msin.datacolumn=${data_column} uvwflag.uvlambdamin=${uvlambdamin} uvwflag.uvlambdamax=${uvlambdamax} msout=${msin}.l${uvlambdamin}to${uvlambdamax} msout.overwrite=True > "${params.data.path}/${params.out.logs}/flag/uvwflag_${msin}_${data_column}_${time}.log" 2>&1
        """
}


process ApplyBEAM {
    label 'sing'

    input:
        val ready
        path ms
        path parset
        val incol
        val outcol

    output:
        val true

    script:
        time = getTime()

        """
        DP3 ${parset} msin=${ms} msin.datacolumn=${incol} msout.datacolumn=${outcol} > "${params.data.path}/${params.out.logs}/calibrate/apply_element_beam_${ms}_${time}.log" 2>&1
        """
}


process ReadMSList {

    input:
        val ready
        val dirname
        val txtname

    output:
        val ms_string, emit: list_str

    exec:
        tlist = file( dirname ).resolve( txtname ).readLines()
        ms_string = tlist.collect {"${it}"}.join(" ")
}
