#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SV = {
   1: [2,6],
   2: [3,7],
   3: [4,8],
   4: [5,9],
   5: [1,9],
   6: [2,10],
   7: [1,8],
   8: [2,9],
   9: [3,10],
  10: [2,3],
  11: [3,4],
  12: [5,6],
  13: [6,7],
  14: [7,8],
  15: [8,9],
  16: [9,10],
  17: [1,4],
  18: [2,5],
  19: [3,6],
  20: [4,7],
  21: [5,8],
  22: [6,9],
  23: [1,3],
  24: [4,6],
  25: [5,7],
  26: [6,8],
  27: [7,9],
  28: [8,10],
  29: [1,6],
  30: [2,7],
  31: [3,8],
  32: [4,9],
}

def shift(register, feedback, output):
    """GPS Shift Register
    
    :param list feedback: which positions to use as feedback (1 indexed)
    :param list output: which positions are output (1 indexed)
    :returns output of shift register:
    
    """
    
    # calculate output
    out = [register[i-1] for i in output]
    if len(out) > 1:
        out = sum(out) % 2
    else:
        out = out[0]
        
    # modulo 2 add feedback
    fb = sum([register[i-1] for i in feedback]) % 2
    
    # shift to the right
    for i in reversed(range(len(register[1:]))):
        register[i+1] = register[i]
        
    # put feedback in position 1
    register[0] = fb
    
    return out

def PRN(n, sv):
    """Build the CA code (PRN) for a given satellite ID
    
    :param int sv: satellite code (1-32)
    :returns list: ca code for chosen satellite
    
    """
    G1 = [1 for i in range(10)]
    G2 = [1 for i in range(10)]

    ca = np.zeros(n)
    for i in range(n):
        g1 = shift(G1, [3,10], [10]) #feedback 3,10, output 10
        g2 = shift(G2, [2,3,6,8,9,10], SV[sv]) #feedback 2,3,6,8,9,10, output 2,6 for sat 1
        ca[i] = ((g1 + g2) % 2)

    ca = ca*2-1
    return ca

f_prn = 10.23e6 / 10  # chipping frequency
chip = 1/f_prn

# init registers

n = 1023
ca_24 = PRN(n,24)
ca_23 = PRN(n,23)

fig_ca, ax_ca = plt.subplots(1,figsize=(10, 4))
t = np.arange(0,n*chip-chip/2, chip)
plt.plot(t,ca_24)
plt.title('ca')
plt.xlabel('C/A chip')
plt.ylim(-1.5,1.5)

ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/chip/1023))
ax_ca.xaxis.set_major_formatter(ticks_x)

fig_auto, ax_auto = plt.subplots(1,figsize=(10, 4))
ca_auto = np.correlate(ca_24, ca_24, 'same')
plt.plot(t,ca_auto)
plt.xlabel('C/A chip')
plt.title('ca auto')
#plt.ylim(-1.5,1.5)

ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/chip/1023))
ax_auto.xaxis.set_major_formatter(ticks_x)

fig_auto_t, ax_auto_t = plt.subplots(1,figsize=(10, 4))
n = 8000
delay = np.zeros(n)
for i in range(0,n):
    delay[i] = i*1023*chip/n

auto_delay =  np.where(np.abs(delay-1023/2*chip) <= chip, 
        1023 - 1023*np.abs(delay-1023/2*chip)/(chip), 
        0)
y = np.interp(delay/chip/1023, t/chip/1023, ca_auto)
plt.plot(delay/chip-0.5*1023,y/np.max(y))
plt.plot(delay/chip-0.5*1023,auto_delay/np.max(auto_delay))
plt.title('Autocorrelation approximation')
plt.xlabel('C/A chip')
plt.xlim(-80,80)

fig_cross, ax_cross = plt.subplots(1,figsize=(10, 4))
ca_cross = np.correlate(ca_24, ca_23, 'same')
plt.plot(t,ca_cross)
plt.title('ca cross')
plt.xlabel('C/A chip')
plt.ylim(-200,200)

ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/chip/1023))
ax_cross.xaxis.set_major_formatter(ticks_x)

plt.show()
