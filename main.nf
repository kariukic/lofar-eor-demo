#!/usr/bin/env nextflow

include {
    makeDirs;
    FlagIntra;
    FilterInter;
    WScleanImage;
    AOFlag;
    Average;
    AOqualityCollect;
    AOqualityCombine;
    WriteMSlist;
    DICalibrate;
    ApplyGains;
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
    mkdir_ch = makeDirs(true) // Always runs

    pre_ch  = shouldRun('PreProcess')     ? PreProcess(mkdir_ch)       : Channel.empty()
    dis_ch  = shouldRun('DISmooth')       ? DISmooth(pre_ch)           : Channel.empty()
    dip_ch  = shouldRun('DIBandpass')     ? DIBandpass(dis_ch)         : Channel.empty()
    avg_ch  = shouldRun('AVG')            ? AVG(dip_ch)                : Channel.empty()
    dd_ch   = shouldRun('DD')             ? DD(avg_ch)                 : Channel.of(true)
    pdd_ch  = shouldRun('PostDD')         ? PostDD(dd_ch)              : Channel.empty()
    ps_ch   = shouldRun('PowerSpectrum')  ? PowerSpectrum(pdd_ch)      : Channel.empty()

    if (shouldRun('FinalImages')) {
        FinalImages(ps_ch)
    }
}


// workflow other {
//     mkdir_ch = makeDirs ( true )
//     pre_ch = PreProcess( mkdir_ch )
//     dis_ch = DISmooth( pre_ch )
//     dip_ch = DIBandpass( dis_ch )
//     avg_ch = AVG( dip_ch )
//     dd_ch = DD( avg_ch )
//     pdd_ch = PostDD( dd_ch )
//     ps_ch = PowerSpectrum ( pdd_ch )
//     FinalImages ( ps_ch )
// }

workflow PreProcess {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath( params.data.ms_files.raw, glob: true, checkIfExists: true, type: 'dir' )
        flag_intra_ch = FlagIntra ( start_ch, mset_ch )

        filter_ch = FilterInter ( flag_intra_ch.collect(), mset_ch )
        
        filtered_mset_ch = filter_ch.collect { "${params.data.path}/" + it.getName().replace( '.FMS.15ch2s.dppp', '.FMS.15ch2s.dppp' ) }

        flag_ch = AOFlag ( true, filtered_mset_ch.flatten(), params.average.lta_to_di.column, params.average.lta_to_di.aoflagger_strategy, 1 )

        averaged_msnames_ch = filter_ch.collect { it.getName().replace( ".FMS.15ch2s.dppp", "_002_3c196.MS" ) }

        all_msets_and_averaged_msnames_ch = filtered_mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        avg_ch = Average ( flag_ch.collect(), all_msets_and_averaged_msnames_ch, params.average.lta_to_di.column, params.average.lta_to_di.timestep, params.average.lta_to_di.freqstep  )

        flagged_filtered_averaged_mset_ch = mset_ch.collect { "${params.data.path}/" + it.getName().replace(  ".FMS.15ch2s.dppp", "_002_3c196.MS"  ) }

        aoq_ch = AOqualityCollect( avg_ch.done_averaging.collect(), flagged_filtered_averaged_mset_ch.flatten(), params.average.lta_to_di.column )

        mslist_ch = WriteMSlist(aoq_ch.collect(), "${params.data.path}/*_002_3c196.MS", params.data.di_mslist)

        AOqualityCombine( mslist_ch.collect(), params.data.di_mslist, "raw_aoq_stats" )

    emit:
        AOqualityCombine.out.done

}


workflow DISmooth {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath("${params.data.path}/${params.data.di_mslist}", checkIfExists: false).splitText().map { file(it.trim()) }

        sols_ch = DICalibrate( start_ch, mset_ch, params.ddecal.di.parset, params.ddecal.di.sourcedb, params.ddecal.di.sols, params.ddecal.di.incol, params.ddecal.di.solint, params.ddecal.di.uvlambdamin, params.ddecal.di.uvlambdamax, params.ddecal.di.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.di.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.di.apply.parset, params.ddecal.di.incol, params.ddecal.di.outcol)

        aoq_ch = AOqualityCollect( true, apply_gains_ch, params.ddecal.di.outcol )
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        wsc_ch = WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.di.outcol )

        AOqualityCombine( wsc_ch.done.collect(), params.data.di_mslist, "di_smooth_aoq_stats" )

    emit:
        AOqualityCombine.out.done

}


workflow DIBandpass {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath("${params.data.path}/${params.data.di_mslist}", checkIfExists: false).splitText().map { file(it.trim()) }

        sols_ch = DIGainCal( start_ch, mset_ch, params.ddecal.bp.parset, params.ddecal.bp.sourcedb, params.ddecal.bp.sols, params.ddecal.bp.incol, params.ddecal.bp.solint, params.ddecal.bp.uvlambdamin, params.ddecal.bp.uvlambdamax, params.ddecal.bp.nchan )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.bp.sols}" }

        mset_and_solutions_ch = mset_ch.merge( all_solutions_ch.flatten() )

        apply_gains_ch = ApplyGains ( sols_ch.collect(), mset_and_solutions_ch, params.ddecal.bp.apply.parset, params.ddecal.bp.incol, params.ddecal.bp.outcol )

        flag_ch = AOFlag ( true, apply_gains_ch, params.ddecal.bp.outcol, params.ddecal.bp.aoflagger_strategy, 0 )

        aoq_ch = AOqualityCollect( flag_ch, apply_gains_ch, params.ddecal.bp.outcol)
    
        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        wsc_ch = WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_subband, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.ddecal.bp.outcol )

        AOqualityCombine( wsc_ch.done.collect(), params.data.di_mslist, "di_bandpass_aoq_stats" )

    emit:
        AOqualityCombine.out.done

}


workflow AVG {

    take:
        start_ch

    main:

        mset_ch = channel.fromPath("${params.data.path}/${params.data.di_mslist}", checkIfExists: false).splitText().map { file(it.trim()) }

        averaged_msnames_ch = mset_ch.collect { it.getName().replace( "_002", "_003" ) }
        
        all_msets_and_averaged_msnames_ch = mset_ch.flatten().merge( averaged_msnames_ch.flatten() )

        avg_ch = Average ( start_ch, all_msets_and_averaged_msnames_ch, params.ddecal.bp.outcol, params.average.ditodd.timestep, params.average.ditodd.freqstep )

        WriteMSlist( avg_ch.done_averaging.collect(), "${params.data.path}/*_003_3c196.MS", params.data.dd_mslist )

    emit:

        WriteMSlist.out

}

workflow DD {
    take:
        start_ch

    main:

        mset_ch = channel.fromPath("${params.data.path}/${params.data.dd_mslist}", checkIfExists: false).splitText().map { file(it.trim()) }

        calibrate_ch = DDCalibrate( start_ch, mset_ch, params.ddecal.dd.parset, params.ddecal.dd.sourcedb, params.ddecal.dd.sols, params.ddecal.dd.incol, params.ddecal.dd.calmode, params.ddecal.dd.solint, params.ddecal.dd.uvlambdamin, params.ddecal.dd.uvlambdamax, params.ddecal.dd.nchan, params.ddecal.dd.usebeam, params.ddecal.dd.beammode, params.ddecal.dd.smoothnessconstraint, params.ddecal.dd.truncateksmoothkernel, params.ddecal.dd.robust_reg, params.ddecal.dd.propagate_sols, params.ddecal.dd.maxiter, params.ddecal.dd.beamproximitylimit, params.ddecal.dd.correctfreqsmearing, params.ddecal.dd.flagstations )

        all_solutions_ch =  mset_ch.collect { it + "/${params.ddecal.dd.sols}" }

        mset_and_sourcedb_ch = mset_ch.flatten().combine( channel.of( params.ddecal.dd.sourcedb ) )

        mset_sourcedb_solutions_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() )

        subtract_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_ch, params.ddecal.dd.subtract.parset, params.ddecal.dd.incol, params.ddecal.dd.outcol )

        aoq_ch = AOqualityCollect( true, subtract_ch, params.ddecal.dd.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName() + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        wsc_ch = WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.ddecal.dd.outcol )

        AOqualityCombine( wsc_ch.done.collect(), params.data.dd_mslist, "dd_smooth_aoq_stats" )

    emit:
        AOqualityCombine.out.done

}


workflow PostDD {
    take:
        start_ch

    main:
        mset_ch = channel.fromPath("${params.data.path}/${params.data.dd_mslist}", checkIfExists: false).splitText().map { file(it.trim()) }

        aoflagger_ch = AOFlag ( start_ch, mset_ch, params.ddecal.dd.outcol, params.postdd.aoflagger_strategy, 1 )

        uvwflag_ch = UVWFlag (aoflagger_ch.collect(), mset_ch, params.ddecal.dd.outcol, params.postdd.uvlambdamin, params.postdd.uvlambdamax)

        beam_ch = ApplyBEAM( uvwflag_ch.collect(), mset_ch, params.postdd.beam.parset, params.ddecal.dd.outcol, params.postdd.beam.outcol )

        aoq_ch = AOqualityCollect( true, beam_ch, params.postdd.beam.outcol )

        im_names_ch =  mset_ch.collect { it.getSimpleName()  + "_" + params.wsclean.imname }

        mset_and_names_ch = mset_ch.merge( im_names_ch.flatten() )

        wsc_ch = WScleanImage ( aoq_ch, mset_and_names_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout_per_timechunk, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit,  params.postdd.beam.outcol )

        AOqualityCombine( wsc_ch.done.collect(), params.data.dd_mslist, "dd_smooth_aoq_stats" )

    emit:
    
        AOqualityCombine.out.done

}

workflow PowerSpectrum {
    take:
        ready

    main:

        ps_dir = "${params.data.path}/${params.out.results}/${params.pspipe.dir}"
        ps_logs_dir = "${params.data.path}/${params.out.logs}/power_spectrum"

        rev_ch = AddRevision( ready, params.pspipe.obsid, params.postdd.beam.outcol, params.data.path, ps_dir, params.pspipe.node, params.pspipe.max_concurrent, params.pspipe.revision, params.pspipe.merge_ms, params.pspipe.aoflag_after_merge_ms, params.pspipe.time_start_index,  params.pspipe.time_end_index )

        RunPSPIPE( ps_dir, rev_ch.toml_file, params.pspipe.obsid, "${params.data.path}/${params.data.dd_mslist}", params.pspipe.merge_ms, params.pspipe.delay_flagger, params.pspipe.vis_flagger, params.pspipe.ml_gpr, params.pspipe.ml_gpr_inj, ps_logs_dir )

    emit:
       RunPSPIPE.out.ready

}


workflow FinalImages {

    take:
        ready

    main:

        mses_sols_ch = ReadMSList( ready, params.data.path, params.data.dd_mslist )

        mses_and_imname_ch = mses_sols_ch.list_str.combine( channel.of( "postdd_corrected" ) )

        WScleanImage ( true, mses_and_imname_ch, params.wsclean.size, params.wsclean.scale, params.wsclean.niter, params.wsclean.pol, params.wsclean.chansout, params.wsclean.minuvl, params.wsclean.maxuvl, params.wsclean.weight, params.wsclean.polfit, params.postdd.beam.outcol )

    emit:

        WScleanImage.out.done

}

// TODO: get rid of all hardcodes
// TODO: Complete modularity/flow combination
// TODO: Run pspipe after every calibration step
// TODO: Make the power spectra for these