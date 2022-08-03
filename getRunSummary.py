#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#   Created By   : Gess Xavier (TAS @ nanoporetech.com)
#   Crated Date  : July 25th 2022
#   version = '1.0'
#-------------------------------------------------------------------------------
""" getRunSummary.py script prints run summary in a tsv format"""
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
import pandas as pd
from colorama import Fore, Back, Style
print(Fore.GREEN + 'getRunSummary.py')
print(Fore.GREEN + 'Author   : Gess Xavier (TAS @ nanoporetech.com)')
print(Fore.GREEN + 'version = 1.0')
print(Fore.GREEN + 'getRunSummary.py script prints run summary in a tsv format')
print(Fore.GREEN + 'Please Note : Only barcodes > 1000 reads will be printed')
print(Style.RESET_ALL)

os.environ['MINKNOW_API_USE_LOCAL_TOKEN'] = '1'
#os.environ['MINKNOW_TRUSTED_CA'] = "/opt/ont/minknow/conf/rpc-certs/ca.crt"

my_parser = argparse.ArgumentParser(description='Run Summary Script')

my_parser.add_argument('--host',
                       '-ho',
                       type=str,
                       help='hostname/IP : localhost')
my_parser.add_argument('--port',
                       '-p',
                       type=str,
                       help='port : 9501')
my_parser.add_argument('--exid',
                       '-e',
                       type=str,
                       help='exid : 20220115-32')
args = my_parser.parse_args()

outfile = open(args.exid+'_summary.tsv', 'w')

try:
    manager = Manager(host=args.host, port=args.port, use_tls=False)
except Exception:
    pass
    manager = Manager(host=args.host, port=9502, use_tls=True)

positions = list(manager.flow_cell_positions())
Headers = ['Machine ID', 'Protocol Run ID', 'Experiment ID', 'Sample ID', 'Flowcell ID', 'Position', 'Run State', 'Start time', 'Estimated bases GB', 'Read Count', 'Read Count in M', 'Basecalled', 'Basecalled Total GB', 'Total Pass GB', 'Total Fail GB', 'Position', 'Barcode', 'Total Reads', 'Pass GB', 'Fail GB', 'Output Path', 'Run time in Days', 'Run time in H']
print('\t'.join(Headers))
outfile.write('\t'.join(Headers)+'\n')

for pos in positions:
    if pos.running:
        connection = pos.connect()
        protocols = connection.protocol.list_protocol_runs()
        for runid in protocols.run_ids:
            run_info = connection.protocol.get_run_info(run_id=runid)
            #print(run_info)
            if args.exid in run_info.user_info.protocol_group_id.value:
            #if exid in run_info.user_info.protocol_group_id.value:
                rrid = run_info.run_id
                machineId = connection.instance.get_machine_id()
                interesting_acquisition_id = run_info.acquisition_run_ids[-1]
                run_info2 = connection.acquisition.get_acquisition_info(run_id=interesting_acquisition_id)
                #print(run_info2)
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
                #print(difftims)
                basecalld = ''
                if run_info2.config_summary.basecalling_enabled is True:
                    basecalld = 'Yes'
                else:
                    basecalld = 'No'
                #print(basecalld)
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



                Res = [machineId.machine_id, str(rrid), str(run_info.user_info.protocol_group_id.value).strip(), str(run_info.user_info.sample_id.value).strip(), run_info.flow_cell.flow_cell_id, str(pos), RunState, str(st), str(round(run_info2.yield_summary.estimated_selected_bases/1000000000, 2)), str(run_info2.yield_summary.read_count), str(round(run_info2.yield_summary.read_count/1000000, 2)), basecalld, basTot, basPass, basFail, '-','-','-','-','-', run_info.output_path, ("%d days, %d hours, %d minutes and %d seconds" % (rd.days, rd.hours, rd.minutes, rd.seconds)), str(round(difftim,2))]
                #print(Res)
                print('\t'.join(Res))
                outfile.write('\t'.join(Res)+'\n')


                # basecalled data extraction
                def barSUM(barcod):
                    stream = connection.statistics.stream_acquisition_output(
                        acquisition_run_id=interesting_acquisition_id,
                        data_selection=minknow_api.statistics_pb2.DataSelection(step = 60 * 60),
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


                    for filter_groups in stream:
                        #print(filter_groups.buckets)
                        for filter_group in filter_groups.buckets:
                            #print(filter_group)
                            #barcodeData = defaultdict(list)
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
                                #print('O', barcodeData[k]['rc'])
                                #print('M', k, max(barcodeData[k]['rc']))
                                maxValDic[k]['rcmax'].append(str(max(barcodeData[k]['rc'])))
                                maxValDic[k]['bpmax'].append(str(max(barcodeData[k]['bp'])))
                                maxValDic[k]['bfmax'].append(str(max(barcodeData[k]['bf'])))



                            for k in maxValDic:
                                if barcod in k:
                                    #print(k)
                                    #print(barcod, '---', barcodeAlias[k])

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
                    #print(res)
                    return res


                x = range(97)
                prevb = '' # step 60 * 60 prints multiple times


                if RunState == 'Complete':
                    for i in x:
                        tabs = '\t'*14
                        res = barSUM('barcode0'+str(i))
                        res2 = barSUM('barcode'+str(i))
                        #print(res)
                        #print(res2)
                        if 'barcode' in res[1]:
                            prevb = res[1].split('\t')[0]
                            if not prevb:
                                TotalReads = res[1].split('\t')[1]
                                if int(TotalReads) > 1000:
                                    GBDataPass = round(int(res[2].strip('\t'))/1000000000 , 2)
                                    GBDataFail = round(int(res[3].strip('\t'))/1000000000 , 2)
                                    #barcAlias = res[4].strip('\t')
                                    fres1 = [tabs, res[0], res[1].split('\t')[0], str(TotalReads), str(GBDataPass), str(GBDataFail)]
                                    print('\t'.join(fres1))
                                    outfile.write('\t'.join(fres1)+'\n')

                        if 'barcode' in res2[1]:
                            prevb = res[1].split('\t')[0]
                            if not prevb:
                                TotalReads = res2[1].split('\t')[1]
                                if int(TotalReads) > 1000:
                                    GBDataPass = round(int(res2[2].strip('\t'))/1000000000 , 2)
                                    GBDataFail = round(int(res2[3].strip('\t'))/1000000000 , 2)
                                    #barcAlias = res2[4].strip('\t')
                                    fres2 = [tabs, res2[0], res2[1].split('\t')[0], str(TotalReads), str(GBDataPass), str(GBDataFail)]
                                    print('\t'.join(fres2))
                                    outfile.write('\t'.join(fres2)+'\n')
