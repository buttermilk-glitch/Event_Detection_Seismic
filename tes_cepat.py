from obspy import UTCDateTime
from obspy.clients.fdsn import Client
c = Client("AUSPASS")
t = UTCDateTime("2021-09-22T00:00:00")
st = c.get_waveforms("S1", "AUKSC", "*", "HHZ", t, t + 3600)
print(st)