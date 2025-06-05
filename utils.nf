#!/usr/bin/env nextflow
// import java.text.SimpleDateFormat

def getTime() {
    def date = new Date()
    def sdf = new java.text.SimpleDateFormat("ddMMyyyy_HHmmss")
    String start_time = sdf.format(date).toString()
    return start_time
}

def shouldRun ( stage ) {
    def workflow_order = [
        'pre-process',
        'di-smooth',
        'di-bandpass',
        'average-for-dd',
        'dd-smooth',
        'post-process',
        'power-spectrum',
        'fullband-image'
    ]
    
    def start_idx = workflow_order.indexOf(params.start_from)
    def stop_idx  = workflow_order.indexOf(params.stop_after)
    def current_idx = workflow_order.indexOf(stage)

    return current_idx >= start_idx && current_idx <= stop_idx
}
