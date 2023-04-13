#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#   Created By    : Gess Xavier (TAS @ nanoporetech.com)
#   Created Date  : Aug 30th 2022
#   Last modified : Apr 13th 2023
#   version = '2.0'
#-------------------------------------------------------------------------------
""" getRunSummary.py script prints run summary in a tsv & xlsx format"""
#-------------------------------------------------------------------------------

import argparse
import sys
import datetime
import dateutil.relativedelta
import os
from minknow_api.manager import Manager
import minknow_api.statistics_pb2
from collections import defaultdict
from collections import Counter
import operator
import ipaddress
import pandas as pd
from colorama import Fore, Back, Style
from tqdm import tqdm
Pversion='v2.0'
print(Fore.GREEN + 'getRunSummary.py')
print(Fore.GREEN + 'Author   : Gess Xavier (TAS @ nanoporetech.com)')
print(Fore.GREEN + 'version = '+Pversion)
print(Fore.GREEN + 'getRunSummary.py script prints run summary in a tsv format')
print(Fore.GREEN + 'Please Note : Only barcodes > 1000 reads will be printed')
print(Fore.BLUE + 'This code uses "tqdm" package, "tqdm" derives from the Arabic word *taqaddum* (تقدّم) which can mean "progress," & "I love you so much" in Spanish (*te quiero demasiado*).')
print(Style.RESET_ALL)

os.environ['MINKNOW_API_USE_LOCAL_TOKEN'] = '1'
#os.environ['MINKNOW_TRUSTED_CA'] = "/opt/ont/minknow/conf/rpc-certs/ca.crt"


def ValidPorts(p):
    port = int(p)
    if port not in (9501, 9502):
        raise argparse.ArgumentTypeError('invalid port! Valid ports : 9501,9502')
    return port

def ValidIP(ip):
    if ip != 'localhost':
        try:
            ips = ipaddress.ip_address(ip)
        except ValueError:
            raise argparse.ArgumentTypeError('Invalid IP address! Try "localhost" or a valid IP')
    #else:
    #    raise argparse.ArgumentTypeError('Invalid IP address! Try localhost or a valid IP')
    return ip


my_parser = argparse.ArgumentParser(description='Run Summary Script')

my_parser.add_argument('--host',
                       '-ho',
                       type=ValidIP,
                       required=True,
                       help='hostname/IP : localhost')
my_parser.add_argument('--port',
                       '-p',
                       type=ValidPorts,
                       required=True,
                       help='port : 9501')
my_parser.add_argument('--exid',
                       '-e',
                       type=str,
                       required=True,
                       help='exid : 20220115-32')
my_parser.add_argument('--version', '-v', action='version', version=Pversion)
args = my_parser.parse_args()

outfile = open(args.exid+'_summary.tsv', 'w')

try:
    manager = Manager(host=args.host, port=args.port)
except Exception:
    pass
    manager = Manager(host=args.host, port=9502, use_tls=True)

positions = list(manager.flow_cell_positions())


HeaderB = ['Machine ID', 'Protocol Run ID', 'Experiment ID', 'Sample ID', 'Flowcell ID', 'Flowcell QC Pass', 'QC PoreCount', 'Position', 'Run State', 'Start time', 'Estimated bases GB', 'Read Count', 'Read Count in M', 'Basecalled', 'Basecalled Total GB', 'Total Pass GB', 'Total Fail GB', 'Position', 'Barcode', 'Alias', 'Total Reads', 'Pass GB', 'Fail GB', 'Output Path', 'Run time in Days', 'Run time in H']
HeaderN = ['Machine ID', 'Protocol Run ID', 'Experiment ID', 'Sample ID', 'Flowcell ID', 'Flowcell QC Pass', 'QC PoreCount','Position', 'Run State', 'Start time', 'Estimated bases GB', 'Read Count', 'Read Count in M', 'Basecalled', 'Basecalled Total GB', 'Total Pass GB', 'Total Fail GB', 'Output Path', 'Run time in Days', 'Run time in H']

chb = 0
chn = 0
pqcResDic = {}
pqcResStatus = ''
pqcResPore = 0
c=len(positions)

cpro = tqdm(total=c, colour="green", position=0, leave=True)
inpro = tqdm(total=0,  position=1, leave=True)
barpro = tqdm(total=0, position=2, leave=True)

for pos in positions:
    cpro.update(1)
    cpro.set_description_str(f'Checking : {pos}')
    if pos.running:
        connection = pos.connect()
        protocols = connection.protocol.list_protocol_runs()
        for runid in protocols.run_ids:
            #print(headcounter)
            run_info = connection.protocol.get_run_info(run_id=runid)
            #print(run_info)
            pqcFID = run_info.pqc_result.flow_cell_id
            pqcPASS = run_info.pqc_result.passed
            pqcPoreCount = run_info.pqc_result.total_pore_count
            pqcResDic[pqcFID] = [pqcFID, pqcPASS, pqcPoreCount]
            #print('PQC-->',pqcFID, pqcPASS, pqcPoreCount)
            if args.exid in run_info.user_info.protocol_group_id.value:
                #for ii in tqdm(range(1), desc = pos):
                PrevHead = ''
            #if exid in run_info.user_info.protocol_group_id.value:
                rrid = run_info.run_id
                machineId = connection.instance.get_machine_id()
                FlowCellID = ''
                if run_info.flow_cell.flow_cell_id:
                    FlowCellID = run_info.flow_cell.flow_cell_id
                else:
                    FlowCellID = run_info.flow_cell.user_specified_flow_cell_id
                if FlowCellID in pqcResDic:
                    #print(pqcResDic[FlowCellID])
                    if pqcResDic[FlowCellID][1] == True:
                        pqcResStatus = 'PASS'
                    else:
                        pqcResStatus = 'FAIL'
                    pqcResPore = pqcResDic[FlowCellID][2]
                else:
                    pqcResStatus = 'QC not found'
                    pqcResPore = 0


                interesting_acquisition_id = ''
                try:
                    interesting_acquisition_id = run_info.acquisition_run_ids[-1]
                    run_info2 = connection.acquisition.get_acquisition_info(run_id=interesting_acquisition_id)
                    #print(run_info2)
                    barcodingOption = run_info2.config_summary.barcoding_enabled
                    if run_info2.state == 1:
                        RunState = 'Active'
                    else:
                        RunState = 'Complete'
                    startTime = str(run_info.start_time.seconds)+'.'+str(run_info.start_time.nanos)
                    endTime = str(run_info.end_time.seconds)+'.'+str(run_info.end_time.nanos)
                    st= datetime.datetime.fromtimestamp(float(startTime))
                    en= datetime.datetime.fromtimestamp(float(endTime))
                    rd = dateutil.relativedelta.relativedelta (en, st)
                    difftim=(en-st).total_seconds()/3600
                    difftims=round((en-st).total_seconds())
                    pos = str(pos).strip(' (running)')
                    proVal = pos+':'+FlowCellID
                    inpro.set_description_str(f'Found : {proVal}')
                    inpro.update(1)
                    basecalld = ''
                    if run_info2.config_summary.basecalling_enabled is True:
                        basecalld = 'Yes'
                    else:
                        basecalld = 'No'
                    basTot = 0
                    basPass = 0
                    basFail = 0
                    if run_info2.yield_summary.basecalled_pass_bases:
                        totBases = run_info2.yield_summary.basecalled_pass_bases + run_info2.yield_summary.basecalled_fail_bases
                        basTot = str(round(totBases/1000000000, 2))
                        basPass = str(round(run_info2.yield_summary.basecalled_pass_bases/1000000000, 2))
                        basFail = str(round(run_info2.yield_summary.basecalled_fail_bases/1000000000, 2))
                    else:
                        basTot = str(basTot)
                        basPass = str(basPass)
                        basFail = str(basFail)
                    Res = [machineId.machine_id, str(rrid), str(run_info.user_info.protocol_group_id.value).strip(), str(run_info.user_info.sample_id.value).strip(), FlowCellID, pqcResStatus, str(pqcResPore), str(pos), RunState, str(st), str(round(run_info2.yield_summary.estimated_selected_bases/1000000000, 2)), str(run_info2.yield_summary.read_count), str(round(run_info2.yield_summary.read_count/1000000, 2)), basecalld, basTot, basPass, basFail, run_info.output_path, ("%d days, %d hours, %d minutes and %d seconds" % (rd.days, rd.hours, rd.minutes, rd.seconds)), str(round(difftim,2))]
                    if barcodingOption == False:
                        chn+=1
                        if chn <=1:
                            #print('\t'.join(HeaderN))
                            outfile.write('\t'.join(HeaderN)+'\n')
                        #print('\t'.join(Res))
                        outfile.write('\t'.join(Res)+'\n')
                    else:
                        chb+=1
                        if chb <=1:
                            #print('\t'.join(HeaderB))
                            outfile.write('\t'.join(HeaderB)+'\n')
                        colsBefore = Res[0:17]
                        colsAfter = Res[17:]

                    # basecalled data extraction
                    def barSUM(barcod):
                        stream = connection.statistics.stream_acquisition_output(
                            acquisition_run_id=interesting_acquisition_id,
                            data_selection=minknow_api.statistics_pb2.DataSelection(step = 4322 * 60),
                            split=minknow_api.statistics_pb2.AcquisitionOutputSplit(
                            alignment_reference=True,
                            barcode_name=True,
                            lamp_barcode_id=True,
                            lamp_target_id=True,
                            ),
                        )
                        barcodeData = defaultdict(list)
                        maxValDic = defaultdict(list)
                        rc = -1
                        rcF = ''
                        bp = -1
                        bpF = ''
                        bf = -1
                        bfF = ''
                        barcodeAlias = {}
                        filtgvals='buckets'
                        filtgval='bucket'

                        for filter_groups in stream:
                            try:
                                for filter_group in filter_groups.buckets:
                                    for grp in filter_group.filtering:
                                        barcode = grp.barcode_name
                                        barcodalias = grp.barcode_alias
                                        barcodeAlias = { barcode: barcodalias}
                                        barcodeData = { barcode: { 'bp': [], 'bf': [], 'rc': [] }  }

                                    for bucket in filter_group.buckets:
                                        bpass = bucket.yield_summary.basecalled_pass_bases
                                        bfail = bucket.yield_summary.basecalled_fail_bases
                                        RcPass = bucket.yield_summary.basecalled_pass_read_count
                                        barcodeData[barcode]['bp'].append(bpass)
                                        barcodeData[barcode]['bf'].append(bfail)
                                        barcodeData[barcode]['rc'].append(RcPass)

                                    for k in barcodeData:
                                        maxValDic = { k: { 'rcmax':[], 'bpmax':[], 'bfmax':[]  } }
                                        maxValDic[k]['rcmax'].append(str(max(barcodeData[k]['rc'])))
                                        maxValDic[k]['bpmax'].append(str(max(barcodeData[k]['bp'])))
                                        maxValDic[k]['bfmax'].append(str(max(barcodeData[k]['bf'])))

                                    for k in maxValDic:
                                        if barcod in k:
                                            if int(maxValDic[k]['rcmax'][0]) > rc:
                                                rc = int(maxValDic[k]['rcmax'][0])
                                                rcF = k+':'+barcodeAlias[k]+'\t'+str(int(maxValDic[k]['rcmax'][0]))

                                            if int(maxValDic[k]['bpmax'][0]) > bp:
                                                bp = int(maxValDic[k]['bpmax'][0])
                                                bpF = '\t'+str(int(maxValDic[k]['bpmax'][0]))

                                            if int(maxValDic[k]['bfmax'][0]) > bf:
                                                bf = int(maxValDic[k]['bfmax'][0])
                                                bfF = '\t'+str(int(maxValDic[k]['bfmax'][0]))
                            except Exception:
                                pass
                                '''
                                I couldn't find another way to make it compatible with minknow_api ver: 4.5.0 & >
                                <= ver 4.5.0 takes buckets/bucket and >= 4.5.0 takes snapshots/snapshot
                                '''
                                for filter_group in filter_groups.snapshots:
                                    for grp in filter_group.filtering:
                                        barcode = grp.barcode_name
                                        barcodalias = grp.barcode_alias
                                        barcodeAlias = { barcode: barcodalias}
                                        barcodeData = { barcode: { 'bp': [], 'bf': [], 'rc': [] }  }

                                    for snapshot in filter_group.snapshots:
                                        bpass = snapshot.yield_summary.basecalled_pass_bases
                                        bfail = snapshot.yield_summary.basecalled_fail_bases
                                        RcPass = snapshot.yield_summary.basecalled_pass_read_count
                                        barcodeData[barcode]['bp'].append(bpass)
                                        barcodeData[barcode]['bf'].append(bfail)
                                        barcodeData[barcode]['rc'].append(RcPass)

                                    for k in barcodeData:
                                        maxValDic = { k: { 'rcmax':[], 'bpmax':[], 'bfmax':[]  } }
                                        maxValDic[k]['rcmax'].append(str(max(barcodeData[k]['rc'])))
                                        maxValDic[k]['bpmax'].append(str(max(barcodeData[k]['bp'])))
                                        maxValDic[k]['bfmax'].append(str(max(barcodeData[k]['bf'])))

                                    for k in maxValDic:
                                        if barcod in k:
                                            if int(maxValDic[k]['rcmax'][0]) > rc:
                                                rc = int(maxValDic[k]['rcmax'][0])
                                                rcF = k+':'+barcodeAlias[k]+'\t'+str(int(maxValDic[k]['rcmax'][0]))

                                            if int(maxValDic[k]['bpmax'][0]) > bp:
                                                bp = int(maxValDic[k]['bpmax'][0])
                                                bpF = '\t'+str(int(maxValDic[k]['bpmax'][0]))

                                            if int(maxValDic[k]['bfmax'][0]) > bf:
                                                bf = int(maxValDic[k]['bfmax'][0])
                                                bfF = '\t'+str(int(maxValDic[k]['bfmax'][0]))

                        res = [pos, rcF, bpF, bfF]
                        return res

                    x = range(97)
                    prevb = '' # step 60 * 60 prints multiple times
                    if RunState == 'Complete':
                        res = []
                        for i in x:
                            res = barSUM('barcode'+str(i).zfill(2))
                            if 'barcode' in res[1]:
                                TotalReads = res[1].split('\t')[1]
                                if int(TotalReads) > 1000:
                                    GBDataPass = round(int(res[2].strip('\t'))/1000000000 , 2)
                                    GBDataFail = round(int(res[3].strip('\t'))/1000000000 , 2)
                                    barc, alias = res[1].split('\t')[0].split(':')
                                    barVal = barc+':'+alias
                                    barpro.set_description_str(f'barcode : {barVal}')
                                    barpro.update(1)
                                    fres1 = ['\t'.join(colsBefore), res[0], barc, alias, str(TotalReads), str(GBDataPass), str(GBDataFail), '\t'.join(colsAfter)]
                                    #print(colsBefore)
                                    #print('\t'.join(fres1))
                                    outfile.write('\t'.join(fres1)+'\n')

                except Exception:
                    pass
outfile.close()
## Read Tsv & Convert to Excel
resFilePath = os.path.join(os.getcwd(), args.exid+'_summary.tsv')
if os.path.exists(resFilePath) and os.stat(resFilePath).st_size != 0:
    resFile = pd.read_csv(resFilePath, sep="\t")
    ir = len(resFile.index)
    ic = len(resFile.columns)
    resFile.to_excel (args.exid+'_summary.xlsx', index=None, header=True)
    tqdm.write(f"Script complete! output stored in")
    tqdm.write(f"tsv format : {resFilePath}")
    tqdm.write(f"xlsx format : {os.path.join(os.getcwd(), args.exid+'_summary.xlsx')}")
else:
    print('Result File Not Created! Did you mistype Experiment ID?')
