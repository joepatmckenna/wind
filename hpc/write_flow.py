import sys
import numpy as np
from datetime import date, datetime, timedelta
import multiprocessing as mp
import timeit
import wind
#print wind.rep.__doc__
#print wind.ev.__doc__
#print wind.sphere_forward_euler.__doc__

HOME_ROOT = '/panfs/storage.local/home-4/jpm12c/'
DATA_ROOT = HOME_ROOT + 'wind/data/'

PI = np.pi
EARTH_RADIUS = 6371.008  # km
MIN_LAT, MAX_LAT = -78.375, 78.375
MIN_LON, MAX_LON = -179.875, 179.875
N_LAT, N_LON = 628, 1440

SECOND = timedelta(seconds=1)
MINUTE = timedelta(minutes=1)
HOUR = timedelta(hours=1)
DAY = timedelta(days=1)
SEC_PER_HR = HOUR.total_seconds()

MIN_DATETIME = datetime(1987, 7, 2)
MAX_DATETIME = datetime(2012, 1, 1)
N_DAYS = (MAX_DATETIME - MIN_DATETIME).days


# return lat/lon files
def lat_lon_data(units='deg'):
    lat = np.linspace(MIN_LAT, MAX_LAT, N_LAT)
    lon = np.linspace(MIN_LON, MAX_LON, N_LON)
    if units == 'rad':
        lat = PI * (lat + 90) / 180
        lon = PI * (lon + 180) / 180
    return lat, lon


def uv_data(date=MIN_DATETIME, units='m_per_s', flat=False):
    filename = DATA_ROOT + date.strftime('%Y/%m/%d/') + '%02i.npz' % (
        date.hour / 6 * 6, )
    with np.load(filename) as data:
        u, v = data['u'], data['v']
    if units == 'km_per_hr':
        u = 3.6 * np.array(u)
        v = 3.6 * np.array(v)
    return u, v


class spline(object):
    def __init__(self,
                 x,
                 y,
                 z,
                 iopt=[0, 0, 0],
                 ider=[-1, 0, -1, 0],
                 z0=None,
                 z1=None,
                 s=0.):

        z = np.concatenate(z)

        nx, tx, ny, ty, c, fp, ier = wind.rep(iopt, ider,
                                              x.copy(),
                                              y.copy(), z.copy(), z0, z1, s)
        if not ier in [0, -1, -2]:
            print 'spline rep ier:', ier

        self.ier = ier
        self.fp = fp
        self.tx = tx[:nx]
        self.ty = ty[:ny]
        self.c = c[:(nx - 4) * (ny - 4)]

    def __call__(self, x, y):

        z, ier = wind.ev(self.tx, self.ty, self.c, x, y)

        if not ier == 0:
            print 'spline ev ier:', ier

        self.ier = ier

        return z

    def ev(self, x, y):
        return self.__call__(x, y)


def sphere_forward_euler((s0, s1, sn, ds, lat0, lon0, u_spl, v_spl)):

    return wind.sphere_forward_euler(
        s0, s1, sn, ds, lat0, lon0, u_spl[0].tx, u_spl[0].ty, u_spl[0].c,
        u_spl[1].tx, u_spl[1].ty, u_spl[1].c, v_spl[0].tx, v_spl[0].ty,
        v_spl[0].c, v_spl[1].tx, v_spl[1].ty, v_spl[1].c)


def integrate(t0,
              tn,
              dt=MINUTE,
              n_lat=10,
              n_lon=10,
              min_lat=MIN_LAT,
              max_lat=MAX_LAT,
              min_lon=MIN_LON,
              max_lon=MAX_LON,
              parallel=True):

    if parallel:
        pool = mp.Pool(processes=mp.cpu_count())

    lat_data, lon_data = lat_lon_data(units='rad')
    lat = np.linspace(np.deg2rad(min_lat), np.deg2rad(max_lat), n_lat)
    lon = np.linspace(np.deg2rad(min_lon), np.deg2rad(max_lon), n_lon)
    lon0, lat0 = np.meshgrid(lon, lat)

    next_data_time = datetime(t0.year, t0.month, t0.day, t0.hour / 6 * 6)
    u_data, v_data = ['', ''], ['', '']
    u_spl, v_spl = ['', ''], ['', '']
    u_data[1], v_data[1] = uv_data(next_data_time, units='km_per_hr')
    u_spl[1] = spline(lat_data, lon_data, u_data[1])
    v_spl[1] = spline(lat_data, lon_data, v_data[1])

    ds = int(dt.total_seconds())
    nt = 1 + int((tn - t0).total_seconds() / ds)
    t = t0 + dt

    flow = np.stack((lat0, lon0), -1)

    t_dur_sec = (tn - t0).total_seconds()
    start_time = timeit.default_timer()

    while t < tn:

        last_data_time = next_data_time
        next_data_time = last_data_time + 6 * HOUR

        u_data[0], v_data[0] = u_data[1], v_data[1]
        u_data[1], v_data[1] = uv_data(next_data_time, units='km_per_hr')

        u_spl[0], v_spl[0] = u_spl[1], v_spl[1]
        u_spl[1] = spline(lat_data, lon_data, u_data[1])
        v_spl[1] = spline(lat_data, lon_data, v_data[1])

        s0 = int((last_data_time - t0).total_seconds())
        s1 = int((t - t0).total_seconds())
        sn = int((min(next_data_time, tn) - t0).total_seconds())
        ns = 1 + int((sn - s1) / ds)

        if parallel:

            args = [(s0, s1, sn, ds, flow[lat_idx, lon_idx, 0],
                     flow[lat_idx, lon_idx, 1], u_spl, v_spl)
                    for lon_idx in range(n_lon) for lat_idx in range(n_lat)]
            data = pool.map_async(sphere_forward_euler, args)
            flow = np.reshape(data.get(), (n_lat, n_lon, 2))

        else:

            for lat_idx in range(n_lat):
                for lon_idx in range(n_lon):
                    flow[lat_idx, lon_idx] = sphere_forward_euler(
                        (s0, s1, sn, ds, flow[lat_idx, lon_idx, 0],
                         flow[lat_idx, lon_idx, 1], u_spl, v_spl))

        t += ns * dt

        print 'integrate: %0.2f%% in %0.2f min' % (
            100 * (t - t0).total_seconds() / t_dur_sec,
            (timeit.default_timer() - start_time) / 60.), datetime.now()

    if parallel:
        pool.close()
        pool.join()

    return flow

if __name__ == '__main__':

    t = sys.argv[1]
    filename = HOME_ROOT + 'wind/vid/data/katrina_lcs/flow/' + str(t) + '.npy'

    t = datetime.strptime(t, '%Y-%m-%d-%H')

    T = 14 * DAY

    flow = integrate(
        t,
        t + T,
        n_lat=768,
        n_lon=768,
        min_lat=MIN_LAT,
        max_lat=MAX_LAT,
        min_lon=MIN_LON,
        max_lon=MAX_LON)

    flow = flow.astype(np.float32)

    np.save(filename, flow)
