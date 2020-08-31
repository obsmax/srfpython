# nice but bad ideaq

# from srfpython.HerrMet.theory import Theory
# from srfpython.HerrMet.parameterizers import Parameterizer
# from srfpython.HerrMet.paramfile import load_paramfile
# from srfpython.HerrMet.datacoders import Datacoder, Datacoder_log, makedatacoder
# import numpy as np
#
# """
#                           a super model array (M)
#   |                          ^    |
#   |                          |    |
#   |                     SuperParameterizer
#   |                          |    |
#   |                          |    v
#   |                         a collection of
#   |                        depth models
#   |                     (ztop, vp, vs, rh)
#   |                            |
#  SuperTheory            Herrmann.HerrmannCaller.disperse
#  (forward problem)             |
#  (Frechet derivatives)         v
#   |                       a collection of
#   |                       dispersion data
#   |                      [(waves, types, modes, freqs, values, (dvalues)), ...]
#   |                          ^    |
#   |                          |    |
#   |                       SuperDatacoder
#   |                          |    |
#   v                          |    v
#                           a super data array (D)
# """
#
#
# # ========== g
# class SuperParameterizer(object):
#
#     def __init__(self, parameterizers):
#         self.parameterizers = parameterizers
#
#     def __len__(self):
#         return len(self.parameterizers)
#
#     def get_Mprior(self):
#         Mprior = []
#         for parameterizer in self.parameterizers:
#             mprior = parameterizer.MMEAN
#             Mprior.append(mprior)
#         return np.hstack(Mprior)
#
#     def split(self, Model):
#         """
#         split a unique Model array into subarrays, one per node
#         :param Model:
#         :return:
#         """
#         j_current = 0
#         models = []
#         for node_number, parameterizer in enumerate(self.parameterizers):
#             nparameters = len(parameterizer.MMEAN)
#             model = Model[j_current : j_current + nparameters]
#             j_current += nparameters
#             models.append(model)
#
#         assert j_current == len(Model)
#         return models
#
#     def inv(self, Model):
#         zvpvsrh_models = []
#         for nnode, model in enumerate(self.split(Model)):
#             ztop, vp, vs, rh = self.parameterizers[nnode].inv(model)
#             zvpvsrh_models.append((ztop, vp, vs, rh))
#         return zvpvsrh_models
#
#     def inv_to_depthmodels(self, Model):
#         dms = []
#         for nnode, model in enumerate(self.split(Model)):
#             dms.append(self.parameterizers[nnode].inv_to_depthmodel(model))
#         return dms
#
#     def inv_to_mod96strings(self, Model):
#         m96strings = []
#         for nnode, model in enumerate(self.split(Model)):
#             m96strings.append(self.parameterizers[nnode].inv_to_mod96string(model))
#         return m96strings
#
#     def show_models(self, ax, Model, Munc=None, offset=3.0, **kwargs):
#         xticks = []
#         xticklabels = []
#         vsticks = np.array([1., 2., 3.])
#         vsticklabels = ['1', '2', '3']
#         if Munc is None:
#             for nnode, dm in enumerate(self.inv_to_depthmodels(Model)):
#                 dm.vs.values += offset * nnode  # offset
#                 dm.vs.show(ax, **kwargs)
#                 xticks = np.concatenate((xticks, vsticks + offset * nnode))
#                 xticklabels = np.concatenate((xticklabels, vsticklabels))
#         else:
#             for nnode, (dmlow, dmhigh) in enumerate(zip(self.inv_to_depthmodels(Model-Munc),
#                                                         self.inv_to_depthmodels(Model+Munc))):
#                 dmlow.vs.values += offset * nnode  # offset
#                 dmhigh.vs.values += offset * nnode  # offset
#                 dmlow.vs.fill_between(ax=ax, other=dmhigh.vs, **kwargs)
#                 xticks = np.concatenate((xticks, vsticks + offset * nnode))
#                 xticklabels = np.concatenate((xticklabels, vsticklabels))
#         ax.set_xticks(xticks)
#         ax.set_xticklabels(xticklabels)
#
#
# class SuperDatacoder(object):
#     def __init__(self, datacoders):
#         self.datacoders = datacoders
#
#     def __len__(self):
#         return len(self.datacoders)
#
#     def get_Dobs(self):
#         Dobs = []
#         for datacoder in self.datacoders:
#             dobs_current, CDinv_diag_current = datacoder.target()
#             Dobs.append(dobs_current)
#         return np.hstack(Dobs)
#
#     def get_Dunc(self, scale_uncertainties=1., add_uncertainty=0.):
#         Dunc = []
#         for datacoder in self.datacoders:
#             dobs_current, CDinv_diag_current = datacoder.target()
#             dunc_current = CDinv_diag_current ** -0.5  # back to uncertainty
#             Dunc.append(scale_uncertainties * dunc_current + add_uncertainty)
#         return np.hstack(Dunc)
#
#     def split(self, Data):
#         i_current = 0
#         datas = []
#         for node_number, datacoder in enumerate(self.datacoders):
#             ndatapoints = len(datacoder.values)
#             data = Data[i_current: i_current + ndatapoints]
#             i_current += ndatapoints
#             datas.append(data)
#
#         assert i_current == len(Data)
#         return datas
#
#     def inv(self, Data):
#         Values = []
#         for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
#             Values.append(datacoder.inv(data))
#         return Values
#
#     def inv_to_laws(self, Data):
#         Laws = []
#         for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
#             Laws.append(datacoder.inv_to_laws(data))
#         return Laws
#
#     def inv_to_surf96strings(self, Data):
#         surf96strings = []
#         for nnode, (datacoder, data) in enumerate(zip(self.datacoders, self.split(Data))):
#             s96 = datacoder.inv_to_surf96string(data)
#             surf96strings.append(s96)
#         return surf96strings
#
#     def show_datas(self, ax, Data, showdvalue=False,
#                    offset=3.,
#                    pticks=[0.5, 1., 2.],
#                    pticklabels=['0.5', '1', '2'],
#                    **kwargs):
#         xticks = []
#         xticklabels = []
#
#         for nnode, laws in enumerate(self.inv_to_laws(Data)):
#
#             xticks = np.concatenate((xticks, np.log(pticks) + offset * nnode))
#             xticklabels = np.concatenate((xticklabels, pticklabels))
#             for law in laws:
#
#                 x = np.log(1. / law.freq) + offset * nnode
#                 y = law.value
#                 if showdvalue:
#                     yinf = law.value * np.exp(-law.dvalue / law.value)
#                     ysup = law.value * np.exp(+law.dvalue / law.value)
#                     ax.fill_between(x, yinf, ysup, **kwargs)
#                 else:
#                     ax.plot(x, y, **kwargs)
#         ax.set_xscale('linear')
#         ax.set_yscale('log')
#         ax.set_xticks(xticks)
#         ax.set_xticklabels(xticklabels)
#
#
# class SuperParameterizerFromStrings(SuperParameterizer):
#     def __init__(self, parameterizer_strings):
#         parameterizers = \
#             [load_paramfile(parameterizer_string, verbose=False)[0]
#              for parameterizer_string in parameterizer_strings]
#
#         SuperParameterizer.__init__(self, parameterizers)
#
#
# class SuperDatacoderFromStrings(SuperDatacoder):
#     def __init__(self, datacoder_strings, which=Datacoder_log):
#         datacoders = \
#             [makedatacoder(datacoder_string, which=which)
#              for datacoder_string in datacoder_strings]
#
#         SuperDatacoder.__init__(self, datacoders)
#
#
# class SuperTheory(object):
#
#     def __init__(self, superparameterizer, superdatacoder):
#         assert len(superparameterizer.list) == len(superdatacoder.list)
#         self.list = [Theory(parameterizer=p, datacoder=d)
#                      for p, d in zip(superparameterizer.list,
#                                      superdatacoder.list)]
#
#     def __call__(self, M, verbose=True):
#         M = M.reshape((nz, ny, nx))
#         theorys = self.theorys
#
#         def job_generator():
#             for iy in range(ny):
#                 for jx in range(nx):  # order matters!!!!
#                     nnode = iy * nx + jx
#                     yield Job(nnode, theorys[iy, jx], M[:, iy, jx])
#
#         def job_handler(nnode, theory, m):
#             data = theory(m=m)
#             return nnode, data
#
#         if verbose:
#             wb = waitbarpipe('g(m)')
#         Data = np.zeros(ny * nx * nper, float)
#         with MapSync(job_handler, job_generator(), **self.mapkwargs) as ma:
#             ib = 0
#             for nnode, (nnode, data), _, _ in ma:
#                 Data[ib: ib + nper] = data
#                 ib += nper
#
#                 if verbose:
#                     wb.refresh(nnode / float(nx * ny))
#         if verbose:
#             wb.close()
#
#         return Data  # warning : Data means encoded data
#
#     def frechet_derivatives(self, M, verbose=True):
#
#         M = M.reshape((nz, ny, nx))
#         theorys = self.theorys
#
#         def job_generator():
#             for iy in range(ny):
#                 for jx in range(nx): # order matters!!!!
#                     nnode = iy * nx + jx
#                     yield Job(nnode, theorys[iy, jx], M[:, iy, jx])
#
#         def job_handler(nnode, theory, m):
#             fd = theory.frechet_derivatives(m=m)
#             return nnode, fd
#
#         if verbose:
#             wb = waitbarpipe('dg/dm(m)')
#
#         # FDs = []
#         rows = []
#         cols = []
#         dats = []
#         with MapSync(job_handler, job_generator(), **self.mapkwargs) as ma:  # order matters
#             for jobid, (nnode, fd), _, _ in ma:
#                 # FDs.append(fd)
#                 iy = nnode // nx  # line number of the node in the surface plane = "latitude number"
#                 ix = nnode % nx  # column number of the node in the surface plane = "longitude number"
#
#                 # _cols, _rows = np.meshgrid(np.arange(nz), np.arange(nper))
#                 # cols += list(nnode * nz + _cols.flat[:])
#                 # rows += list(nnode * nper + _rows.flat[:])
#
#                 # iz = depth number
#                 iz, ida = np.meshgrid(np.arange(nz), np.arange(nper))
#                 # convert indexs into positions in the model (_cols) ad data (_rows) spaces
#                 _cols = iz.flat[:] * (nx * ny) + iy * nx + ix
#                 _rows = list(nnode * nper + ida.flat[:])
#                 rows += list(_rows)
#                 cols += list(_cols)
#                 dats += list(fd.flat[:])
#
#                 if verbose:
#                     wb.refresh(nnode / float(nx * ny))
#
#         if verbose:
#             wb.close()
#
#         G = sp.csc_matrix((dats, (rows, cols)), shape=(nper*ny*nx, nz*ny*nx), dtype=float)
#         return G
#
