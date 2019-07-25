# Here we test the basic functionality of reading an ensemble from
# a TRDI ADCP binary output file (PD0 format)
#
# these tests are designed to follow the example here
# http://doc.pytest.org/en/latest/fixture.html
import pytest


@pytest.fixture
def trdipd0_analyzepd0file():
    import TRDIstuff.TRDIpd0tonetcdf as TRDIpd0
    # this is a test data file, very simple ADCP file
    pd0file = '9991wh000.000'
    # this command will display basics about the raw binary data
    maxens, enslen, ensdata, startofdata = TRDIpd0.analyzepd0file(pd0file)
    return maxens, enslen, ensdata, startofdata


def test__maxens(trdipd0_analyzepd0file):
    import datetime
    maxens, enslen, ensdata, startofdata = trdipd0_analyzepd0file
    assert maxens == 500
    assert enslen == 468
    assert startofdata == 0
    assert ensdata['Header'] == {'headerID': 127, 'sourceID': 127,
                                 'nbytesperens': 466, 'ndatatypes': 6,
                                 'offsets': [18, 71, 136, 266, 332, 398]}
    assert ensdata['FLeader'] == {'CPU_Version': '16.202',
                                  'System_Configuration_LSB': '11001010',
                                  'System_Frequency': 300, 'Beam_Pattern': 'Convex',
                                  'Sensor_Configuration': 1,
                                  'Transducer_Head_Is_Attached': 'Yes',
                                  'Orientation': 'Up-facing beams',
                                  'System_Configuration_MSB': '01000001',
                                  'Beam_Angle': 20, 'Beam_Configuration': '4-bm janus',
                                  'Simulated_Data': 0, 'Lag_Length': 53,
                                  'Number_of_Beams': 4, 'Number_of_Cells': 16,
                                  'Pings_Per_Ensemble': 300, 'Depth_Cell_Length_cm': 200,
                                  'Blank_after_Transmit_cm': 176,
                                  'Signal_Processing_Mode': 1,
                                  'Low_Corr_Threshold': 64, 'No._Code_Reps': 5,
                                  'PGd_Minimum': 0, 'Error_Velocity_Threshold': 2000,
                                  'Time_Between_Ping Groups': '000:01:00',
                                  'Coord_Transform_LSB': '00000111',
                                  'Coord_Transform': 'BEAM', 'Tilts_Used': 'Yes',
                                  '3-Beam_Solution_Used': 'Yes',
                                  'Bin_Mapping_Used': 'Yes',
                                  'Heading_Alignment_Hundredths_of_Deg': 0,
                                  'Heading_Bias_Hundredths_of_Deg': 0,
                                  'Sensor_Source_Byte': '01111111',
                                  'Calculate_EC_from_ED_ES_and_ET': 'Yes',
                                  'Uses_ED_from_depth_sensor': 'Yes',
                                  'Uses_EH_from_transducer_heading_sensor': 'Yes',
                                  'Uses_EP_from_transducer_pitch_sensor': 'Yes',
                                  'Uses_ER_from_transducer_roll_sensor': 'Yes',
                                  'Uses_ES_from_conductivity_sensor': 'Yes',
                                  'Uses_ET_from_transducer_temperature_sensor': 'Yes',
                                  'Sensor_Avail_Byte': '00011101',
                                  'Speed_of_sound_sensor_available': 'No',
                                  'Depth_sensor_available': 'No',
                                  'Heading_sensor_available': 'Yes',
                                  'Pitch_sensor_available': 'Yes',
                                  'Roll_sensor_available': 'Yes',
                                  'Conductivity_sensor_available': 'No',
                                  'Temperature_sensor_available': 'Yes',
                                  'Bin_1_distance_cm': 422, 'Xmit_pulse_length_cm': 244,
                                  'Ref_Lyr_Avg_Starting_cell': 1,
                                  'Ref_Lyr_Avg_Ending_cell': 5,
                                  'False_Target_Threshold': 50,
                                  'Transmit_lag_distance_cm': 49,
                                  'CPU_Board_Serial_Number': '8800099e1a9',
                                  'System_Bandwidth': 0, 'System_Power': 255,
                                  'Base_Frequency_Index': 128}
    # note that these needed to be tested separately because of the datetime conversion
    # in this part of the data
    assert ensdata['VLeader'] == {'Ensemble_Number': 200, 'Year': 2004, 'Month': 5,
                                  'Day': 19, 'Hour': 2, 'Minute': 42, 'Second': 30,
                                  'Hundredths': 0, 'Ensemble_#_MSB': 0,
                                  'timestr': '2004:05:19 02:42:30.000',
                                  'dtobj': datetime.datetime(2004, 5, 19, 2, 42, 30),
                                  'julian_day_from_as_datetime_object': 2453144.612847222,
                                  'julian_day_from_julian': 2453421.112847222,
                                  'EPIC_time': 2453421, 'EPIC_time2': 9749999,
                                  'BIT_Result_Byte_13': '00000000',
                                  'Demod_1_error_bit': 0, 'Demod_0_error_bit': 0,
                                  'Timing_Card_error_bit': 0, 'Speed_of_Sound': 1525,
                                  'Depth_of_Transducer': 320,
                                  'Heading, Pitch, Roll units': 'hundredths_of_a_degree',
                                  'Heading': 25701, 'Pitch': -355, 'Roll': 154,
                                  'Salinity': 35, 'Temperature': 2107, 'MPT_minutes': 0,
                                  'MPT_seconds': 0, 'MPT_hundredths': 5,
                                  'H/Hdg_Std_Dev': 0, 'P/Pitch_Std_Dev': 0,
                                  'R/Roll_Std_Dev': 0, 'Xmit_Current': 53,
                                  'Xmit_Voltage': 144, 'Ambient_Temp': 87,
                                  'Pressure_(+)': 255, 'Pressure_(-)': 255,
                                  'Attitude_Temp': 87, 'Attitude': 129,
                                  'Contamination_Sensor': 159,
                                  'Error_Status_Word_Low_16_bits_LSB': '00000000',
                                  'Bus_Error_exception': 0, 'Address_Error_exception': 0,
                                  'Illegal_Instruction_exception': 0,
                                  'Zero_Divide_exception': 0, 'Emulator_exception': 0,
                                  'Unassigned_exception': 0,
                                  'Watchdog_restart_occurred': 0,
                                  'Battery_Saver_power': 0,
                                  'Error_Status_Word_Low_16_bits_MSB': '10000001',
                                  'Pinging': 1, 'Cold_Wakeup_occurred': 0,
                                  'Unknown_Wakeup_occurred': 1,
                                  'Error_Status_Word_High_16_bits_LSB': '00000000',
                                  'Clock_Read_error_occurred': 0, 'Unexpected_alarm': 0,
                                  'Clock_jump_forward': 0, 'Clock_jump_backward': 0,
                                  'Error_Status_Word_High_16_bits_MSB': '00000000',
                                  'Power_Fail_(Unrecorded)': 0,
                                  'Spurious_level_4_intr_(DSP)': 0,
                                  'Spurious_level_5_intr_(UART)': 0,
                                  'Spurious_level_6_intr_(CLOCK)': 0,
                                  'Level_7_interrupt_occurred': 0,
                                  'Pressure_deca-pascals': 0,
                                  'Pressure_variance_deca-pascals': 0,
                                  'RTC_Century': 20, 'RTC_Year': 4, 'RTC_Month': 5,
                                  'RTC_Day': 19, 'RTC_Hour': 2, 'RTC_Minute': 42,
                                  'RTC_Second': 30, 'RTC_Hundredths': 0}
    # note the use of (expression to evaluate).any() to test the array data
    # http://justinbois.github.io/bootcamp/2016/lessons/l19_numpy_arrays.html
    assert (ensdata['VData'] == [
        [-32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
         -32768, -32768, -32768, -32768, -32768],
        [-149, -185, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
         -32768, -32768, -32768, -32768],
        [-32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
         -32768, -32768, -32768, -32768, -32768],
        [-32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
         -32768, -32768, -32768, -32768, -32768]]).all()
    assert (ensdata['CData'] == [[14, 12, 11, 12, 12, 11, 12, 12, 11, 11, 12, 11, 11, 11, 12, 11],
                                 [72, 49, 24, 12, 11, 12, 11, 12, 12, 11, 12, 11, 12, 11, 12, 12],
                                 [23, 12, 11, 12, 11, 11, 11, 12, 11, 11, 11, 12, 11, 11, 11, 12],
                                 [15, 12, 12, 12, 12, 11, 12, 12, 13, 12, 12, 12, 12, 12, 12, 13]]).all()
    assert (ensdata['IData'] == [[40, 40, 40, 39, 40, 39, 39, 39, 40, 39, 39, 39, 39, 39, 39, 39],
                                 [57, 46, 43, 43, 43, 43, 43, 43, 43, 43, 42, 43, 43, 43, 43, 43],
                                 [40, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39],
                                 [41, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39]]).all()
    assert (ensdata['GData'] == [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [92, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]).all()
