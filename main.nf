#!/usr/bin/env nextflow

include {
    setupDirs;
    FlagIntra;
    FilterInter;
    WScleanImage;
    AOFlag;
    AverageData;
    AOqualityCollect;
    AOqualityCombine;
    GetData;
    DICalibrate;
    ApplyGains;
    H5ParmCollect;
    DIGainCal;
    DDCalibrate;
    SubtractSources;
    UVWFlag;
    ApplyBEAM;
    ReadMSList;
} from './processes.nf'

include {
    AddRevision;
    RunPSPIPE
} from "./makeps.nf"

include {
    shouldRun
} from './utils.nf'


workflow {
    dir_ch = setupDirs(true) // Always runs
    pre_ch  = shouldRun('pre-process')      ? PreProcess(dir_ch)    : true //Channel.of(true)
    dis_ch  = shouldRun('di-smooth')        ? DISmooth(pre_ch)      : true //Channel.of(true)
    dip_ch  = shouldRun('di-bandpass')      ? DIBandpass(dis_ch)    : true //Channel.of(true)
    avg_ch  = shouldRun('average-for-dd')   ? AverageForDD(dip_ch)  : true //Channel.of(true)
    dds_ch   = shouldRun('dd-smooth')       ? DDSmooth(avg_ch)      : true //Channel.of(true)
    pdd_ch  = shouldRun('post-process')     ? PostProcess(dds_ch)   : true //Channel.of(true)
    ps_ch   = shouldRun('power-spectrum')   ? PowerSpectrum(pdd_ch) : true //Channel.of(true)
            if (shouldRun('fullband-image')) {
                Image(ps_ch)
            }
    }


// workflow {
//     dir_ch = setupDirs(true)
//     // pre_ch = preprocess( dir_ch )
//     // dis_ch = DISmooth( pre_ch )
//     // dip_ch = DIBandpass( dis_ch )
//     // avg_ch = AverageForDD( dip_ch )
//     // dd_ch = DDSmooth( dir_ch ) //avg_ch )
//     // pdd_ch = PostProcess( dd_ch )
//     // ps_ch = PowerSpectrum ( pdd_ch )
//     // FinalImages ( ps_ch )
// }


workflow PreProcess {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( params.data.ms_files.raw, glob: true, checkIfExists: true, type: 'dir' )
        flag_intra_ch = FlagIntra ( start_ch, mset_ch )

        filter_ch = FilterInter ( flag_intra_ch.collect(), mset_ch )
        
        filtered_mset_ch = filter_ch.collect { "${params.data.path}/" + it.getName() }  //.replace( '.FMS.15ch2s.dppp', '.FMS.15ch2s.dppp' ) }

        flag_ch = AOFlag ( true, filtered_mset_ch.flatten(), params.average.lta_to_di.column, params.average.lta_to_di.aoflagger_strategy, 1 )

        averaged_msnames_ch = filter_ch.collect { it.getName() + '.flagged.di_averaged' } //.replace( ".FMS.15ch2s.dppp", "_002_3c196.MS" ) }

        all_msets_and_averaged_msnames_ch = filtered_mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        avg_ch = AverageData ( flag_ch.collect(), all_msets_and_averaged_msnames_ch, params.average.lta_to_di.column, params.average.lta_to_di.timestep, params.average.lta_to_di.freqstep  )

        flagged_filtered_averaged_mset_ch = mset_ch.collect { "${params.data.path}/" + it.getName() + '.noInter.flagged.di_averaged'  } //.replace( ".FMS.15ch2s.dppp", "_002_3c196.MS"  ) }

        aoq_ch = AOqualityCollect( avg_ch.done_averaging.collect(), flagged_filtered_averaged_mset_ch.flatten(), params.average.lta_to_di.column )

        mslist_ch = GetData(aoq_ch.collect(), "${params.data.path}/*noInter.flagged.di_averaged", params.data.di_mslist)

        AOqualityCombine( mslist_ch.collect(), "${params.data.path}/${params.data.di_mslist}", "raw_aoq_stats" )

    emit:
        AOqualityCombine.out.done

}


workflow DISmooth {

    take:
        start_ch

    main:

        // mset_ch = GetData(start_ch, "${params.data.path}/*noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file( it.trim( ) ) }

        def mslist_file = file(params.data.di_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file(it.trim()) }
        }

        sols_ch = DICalibrate( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint, params.ddecal.di.uvlambdamin, params.ddecal.di.uvlambdamax, params.ddecal.di.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.di.sols}" }

        sols_collect_ch = H5ParmCollect( sols_ch.collect(), all_solutions_ch.collect().map{ it.join(' ') }, "di_smooth_solutions")

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_collect_ch.done, mset_and_solutions_ch, params.ddecal.di.apply.parset, params.ddecal.di.incol, params.ddecal.di.outcol)

        aoq_ch = AOqualityCollect( true, apply_gains_ch, params.ddecal.di.outcol )

        aoq_comb_ch = AOqualityCombine( aoq_ch.collect(), "${params.data.path}/${params.data.di_mslist}", "di_smooth_aoq_stats" )

        ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}_di"
        ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum/di"

        rev_ch = AddRevision( aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.di.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index  )

        ps_ch = RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.di_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( ps_ch.ready, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.di.outcol )

    emit:
        WScleanImage.out.done.collect()

}


workflow DIBandpass {

    take:
        start_ch

    main:

        // mset_ch = GetData(start_ch, "${params.data.path}/*noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file( it.trim( ) ) }

        def mslist_file = file(params.data.di_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file(it.trim()) }
        }

        sols_ch = DIGainCal( start_ch, mset_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint, params.ddecal.bp.uvlambdamin, params.ddecal.bp.uvlambdamax, params.ddecal.bp.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        sols_collect_ch = H5ParmCollect( sols_ch.collect(), all_solutions_ch.collect().map{ it.join(' ') }, "di_bandpass_solutions")

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_collect_ch.done, mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol )

        flag_ch = AOFlag ( true, apply_gains_ch, params.ddecal.bp.outcol, params.ddecal.bp.aoflagger_strategy, 0 )

        aoq_ch = AOqualityCollect( flag_ch, apply_gains_ch, params.ddecal.bp.outcol)

        aoq_comb_ch = AOqualityCombine( aoq_ch.collect(), "${params.data.path}/${params.data.di_mslist}", "di_bandpass_aoq_stats" )

        ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}_bp"
        ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum/bp"

        rev_ch = AddRevision( aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.bp.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index  )

        ps_ch = RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.di_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( ps_ch.ready, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_subband, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.bp.outcol )

    emit:
        WScleanImage.out.done.collect()

}


workflow AverageForDD {

    take:
        start_ch

    main:

        // mset_ch = channel.fromPath("${params.data.path}/${params.data.di_mslist}", checkIfExists: false).splitText().map { file( it.trim( ) ) }
        // mset_ch = GetData(start_ch, "${params.data.path}/*noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file( it.trim( ) ) }

        def mslist_file = file(params.data.di_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.noInter.flagged.di_averaged", params.data.di_mslist).splitText().map { file(it.trim()) }
        }

        averaged_msnames_ch = mset_ch.collect { it.getName() + '.dd_averaged' } //.replace( ".di_averaged", ".dd_averaged" ) }
        
        all_msets_and_averaged_msnames_ch = mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        AverageData ( start_ch, all_msets_and_averaged_msnames_ch, params.ddecal.bp.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep )

    emit:

        AverageData.out.done_averaging.collect()

}


workflow DDSmooth {
    take:
        start_ch

    main:

        // mset_ch = channel.fromPath("${params.data.path}/${params.data.dd_mslist}", checkIfExists: false).splitText().map { file( it.trim( ) ) }
        // mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged", params.data.dd_mslist).splitText().map { file( it.trim( ) ) }

        def mslist_file = file(params.data.dd_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged", params.data.dd_mslist).splitText().map { file(it.trim()) }
        }

        calibrate_ch = DDCalibrate( start_ch, mset_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.calmode, params.ddecal.dd.solint, params.ddecal.dd.uvlambdamin, params.ddecal.dd.uvlambdamax, params.ddecal.dd.nchan, params.ddecal.dd.usebeam, params.ddecal.dd.beammode, params.ddecal.dd.smoothnessconstraint, params.ddecal.dd.truncateksmoothkernel, params.ddecal.dd.robust_reg, params.ddecal.dd.propagate_sols, params.ddecal.dd.maxiter, params.ddecal.dd.beamproximitylimit, params.ddecal.dd.correctfreqsmearing, params.ddecal.dd.flagstations )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.dd.sols}" }

        sols_collect_ch = H5ParmCollect( calibrate_ch.collect(), all_solutions_ch.collect().map{ it.join(' ') }, "dd_smooth_solutions")

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.dd.sourcedb ) )

        mset_sourcedb_solutions_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() )

        subtract_ch = SubtractSources ( sols_collect_ch.done, mset_sourcedb_solutions_ch, params.ddecal.dd.subtract.parset, params.ddecal.dd.incol, params.ddecal.dd.outcol )

        aoflagger_ch = AOFlag ( true, subtract_ch, params.ddecal.dd.outcol, params.postdd.aoflagger_strategy, 1 )

        aoq_ch = AOqualityCollect( aoflagger_ch, mset_ch, params.ddecal.dd.outcol )

        aoq_comb_ch = AOqualityCombine( aoq_ch.collect(), "${params.data.path}/${params.data.dd_mslist}", "dd_smooth_aoq_stats" )

        ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}_dd"
        ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum/dd"

        rev_ch = AddRevision( aoq_comb_ch.done, params.pspipe.obsid, params.ddecal.dd.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index  )

        ps_ch = RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.dd_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( ps_ch.ready, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.dd.outcol )

    emit:
        WScleanImage.out.done.collect()

}


workflow PostProcess {
    take:
        start_ch

    main:
        
        // mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged", params.data.dd_mslist).splitText().map { file( it.trim( ) ) }


        def mslist_file = file(params.data.dd_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged", params.data.dd_mslist).splitText().map { file(it.trim()) }
        }

        beam_ch = ApplyBEAM( start_ch, mset_ch, params.postdd.beam.parset, params.ddecal.dd.outcol, params.postdd.beam.outcol )

        UVWFlag ( beam_ch.collect(), mset_ch, params.postdd.beam.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax )

    emit:
        UVWFlag.out.done.collect()

}


workflow PowerSpectrum {
    take:
        start_ch

    main:

        def pdd_ext = "l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}"
        // mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged.${pdd_ext}", params.data.pdd_mslist).splitText().map { file( it.trim( ) ) }
        def mslist_file = file(params.data.pdd_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged.${pdd_ext}", params.data.pdd_mslist).splitText().map { file(it.trim()) }
        }

        aoq_ch = AOqualityCollect( start_ch, mset_ch, params.pspipe.incol )

        aoq_comb_ch = AOqualityCombine( aoq_ch.collect(), "${params.data.path}/${params.data.pdd_mslist}", "pdd_smooth_aoq_stats" )

        ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}_pdd"
        ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum/pdd"

        rev_ch = AddRevision( aoq_comb_ch.done, params.pspipe.obsid, params.pspipe.incol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index  )

        ps_ch = RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.pdd_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

        im_names_ch =  mset_ch.map { it.getSimpleName() +  "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        WScleanImage ( ps_ch.ready, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.pspipe.incol )

    emit:
        WScleanImage.out.done.collect()

}


workflow Image {

    take:
        start_ch

    main:

        // mses_sols_ch = ReadMSList( ready, params.data.path, params.data.pdd_mslist )
        // mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "postdd_corrected_fullband" ) )

        def pdd_ext = "l${params.postdd.uvlambdamin}to${params.postdd.uvlambdamax}"
        def mslist_file = file(params.data.pdd_mslist)
        if (mslist_file.exists()) {
            mset_ch = Channel.value(mslist_file).splitText().map { file(it.trim()) }
        } else {
            mset_ch = GetData(start_ch, "${params.data.path}/*.dd_averaged.${pdd_ext}", params.data.pdd_mslist).splitText().map { file(it.trim()) }
        }

        mses_and_imname_ch = mset_ch.collect().map{ it.join(' ') }.combine( channel.of( "postdd_corrected_fullband" ) )

        WScleanImage ( start_ch, mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.pspipe.incol )

    emit:

        WScleanImage.out.done

}


// workflow PowerSpectrum2 {
//     take:
//         ready

//     main:

//         ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}"
//         ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum"

//         rev_ch = AddRevision( ready, params.pspipe.obsid, params.postdd.beam.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index )

//         RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.dd_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

//     emit:
//        RunPSPIPE.out.ready

// }

// TODO: get rid of all hardcodes
// TODO: Complete modularity/flow combination
// TODO: Run pspipe after every calibration step
// TODO: Make the power spectra for these