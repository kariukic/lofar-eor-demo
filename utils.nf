#!/usr/bin/env nextflow
import java.text.SimpleDateFormat

def getTime() {
    def date = new Date()
    def sdf = new SimpleDateFormat("dd_MM_yyyy_HH_mm_ss")
    String start_time = sdf.format(date).toString() //workflow.start.format("yyyy_MM_dd_HH_mm_ss")
}

def shouldRun ( stage ) {
    def workflow_order = [
        'PreProcess',
        'DISmooth',
        'DIBandpass',
        'AVG',
        'DD',
        'PostDD',
        'PowerSpectrum',
        'FinalImages'
    ]
    
    def start_idx = workflow_order.indexOf(params.start_from)
    def stop_idx  = workflow_order.indexOf(params.stop_after)
    def current_idx = workflow_order.indexOf(stage)

    return current_idx >= start_idx && current_idx <= stop_idx
}
