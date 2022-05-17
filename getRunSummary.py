import argparse
import sys
import datetime
import dateutil.relativedelta
import os
from minknow_api.manager import Manager

#os.environ['MINKNOW_API_USE_LOCAL_TOKEN'] = '1'
os.environ['MINKNOW_TRUSTED_CA'] = "/opt/ont/minknow/conf/rpc-certs/ca.crt"

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

manager = Manager(host=args.host, port=args.port)
positions = list(manager.flow_cell_positions())
Headers = ['Machine ID', 'Protocol run ID', 'Experiment ID', 'Sample ID', 'Flowcell ID', 'Position', 'Run State', 'Start time', 'Estimated bases GB', 'Read Count', 'Read Count in M', 'Output Path', 'Run time in Days', 'Run time in H']
print('\t'.join(Headers))

for pos in positions:
    if pos.running:
        connection = pos.connect()
        protocols = connection.protocol.list_protocol_runs()
        for runid in protocols.run_ids:
            run_info = connection.protocol.get_run_info(run_id=runid)
            if args.exid in run_info.user_info.protocol_group_id.value:
                rrid = run_info.run_id
                machineId = connection.instance.get_machine_id()
                interesting_acquisition_id = run_info.acquisition_run_ids[-1]
                run_info2 = connection.acquisition.get_acquisition_info(run_id=interesting_acquisition_id)
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
                pos = str(pos).strip(' (running)')
                Res = [machineId.machine_id, str(rrid), str(run_info.user_info.protocol_group_id.value).strip(), str(run_info.user_info.sample_id.value).strip(), run_info.flow_cell.flow_cell_id, str(pos), RunState, str(st), str(round(run_info2.yield_summary.estimated_selected_bases/1000000000, 2)), str(run_info2.yield_summary.read_count), str(round(run_info2.yield_summary.read_count/1000000, 2)), run_info.output_path, ("%d days, %d hours, %d minutes and %d seconds" % (rd.days, rd.hours, rd.minutes, rd.seconds)), str(round(difftim,2))]
                print('\t'.join(Res))
