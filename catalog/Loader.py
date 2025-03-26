import numpy as np
from lupa import LuaRuntime
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

class halo(object):
  def __init__(self, path):
    self.Lbox, self.redshift, self.nc, self.omega_m, self.h, self.sigma8, self.fof, self.snp = load_params(path)

    self.position, self.velocity, self.mass = load_data(self)

def load_params(path):
  lua = LuaRuntime(unpack_returned_tuples=True)
  with open(path, "r", encoding="utf-8") as f:
    lua_code = f.read()

  # Lua ファイルを実行
  lua.execute(lua_code)

  # 変数を取得
  nc = lua.eval("nc")
  boxsize = lua.eval("boxsize")
  omega_m = lua.eval("omega_m")
  sigma8 = lua.eval("sigma8")
  h = lua.eval("h")
  output_redshifts = lua.eval("output_redshifts")
  redshift = output_redshifts[0]
  read_powerspectrum = lua.eval("read_powerspectrum")
  write_snapshot = lua.eval("write_snapshot")
  fof = lua.eval("fof")
  snp = lua.eval("snp")

  print(f"nc: {nc}, boxsize: {boxsize}")
  print(f"omega_m: {omega_m}, h: {h}")
  print(f"output_redshifts: {output_redshifts}")
  print(f"read_powerspectrum: {read_powerspectrum}")
  print(f"write_snapshot: {write_snapshot}")
  return boxsize, redshift, nc, omega_m, h, sigma8, fof, snp

def load_data(self):
    m_particle = particle_mass(self)

    

  
  return position, velocity, mass

def particle_mass(self):
    Omega_m = self.omega_m  # 現在の物質密度パラメータ
    box_size = (self.Lbox * u.Mpc)**3  # 500 Mpc/h のシミュレーションボックス (comoving)
    N_particles = self.nc**3  # 総粒子数 (1024^3 の格子)
    H0 = self.h*100

    cosmo = FlatLambdaCDM(H0=H0, Om0=Omega_m)
    rho_crit = cosmo.critical_density0 # g/cm^3

    # 物質の平均密度
    rho_m = Omega_m * rho_crit  # g/cm^3

    # ボックス内の総物質質量
    M_total = rho_m * box_size  # g

    # 粒子質量
    m_particle = (M_total / N_particles).to(u.Msun)  # 太陽質量単位に変換
    return mparticle
    