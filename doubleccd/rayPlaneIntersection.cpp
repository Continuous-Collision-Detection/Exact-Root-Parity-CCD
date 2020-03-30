#include "implicit_point.h"
#include "rayPlaneIntersection.h"
#pragma intrinsic(fabs)

int deltaVolume_filtered(double px, double py, double pz, double dx, double dy, double dz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   double qx_px = qx - px;
   double qy_py = qy - py;
   double rx_px = rx - px;
   double ry_py = ry - py;
   double rz_pz = rz - pz;
   double qz_pz = qz - pz;
   double sx_px = sx - px;
   double sy_py = sy - py;
   double sz_pz = sz - pz;
   double tmp_a = qx_px * ry_py;
   double tmp_b = qy_py * rx_px;
   double m01 = tmp_a - tmp_b;
   double tmq_a = qx_px * rz_pz;
   double tmq_b = qz_pz * rx_px;
   double m02 = tmq_a - tmq_b;
   double tmr_a = qy_py * rz_pz;
   double tmr_b = qz_pz * ry_py;
   double m12 = tmr_a - tmr_b;
   double mt1 = m01 * sz_pz;
   double mt2 = m02 * sy_py;
   double mt3 = m12 * sx_px;
   double mtt = mt1 - mt2;
   double m012 = mtt + mt3;
   double qx_px1 = qx_px - dx;
   double qy_py1 = qy_py - dy;
   double rx_px1 = rx_px - dx;
   double ry_py1 = ry_py - dy;
   double rz_pz1 = rz_pz - dz;
   double qz_pz1 = qz_pz - dz;
   double sx_px1 = sx_px - dx;
   double sy_py1 = sy_py - dy;
   double sz_pz1 = sz_pz - dz;
   double tmp_a1 = qx_px1 * ry_py1;
   double tmp_b1 = qy_py1 * rx_px1;
   double m011 = tmp_a1 - tmp_b1;
   double tmq_a1 = qx_px1 * rz_pz1;
   double tmq_b1 = qz_pz1 * rx_px1;
   double m021 = tmq_a1 - tmq_b1;
   double tmr_a1 = qy_py1 * rz_pz1;
   double tmr_b1 = qz_pz1 * ry_py1;
   double m121 = tmr_a1 - tmr_b1;
   double mt11 = m011 * sz_pz1;
   double mt21 = m021 * sy_py1;
   double mt31 = m121 * sx_px1;
   double mtt1 = mt11 - mt21;
   double m0121 = mtt1 + mt31;
   double res = m012 - m0121;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(dx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(dz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qy_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ry_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rz_pz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qz_pz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sy_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sz_pz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= 5.551115123125786e-14;
   if (res > epsilon) return IP_Sign::POSITIVE;
   if (-res > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int deltaVolume_interval(interval_number px, interval_number py, interval_number pz, interval_number dx, interval_number dy, interval_number dz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz)
{
   setFPUModeToRoundUP();
   interval_number qx_px(qx - px);
   interval_number qy_py(qy - py);
   interval_number rx_px(rx - px);
   interval_number ry_py(ry - py);
   interval_number rz_pz(rz - pz);
   interval_number qz_pz(qz - pz);
   interval_number sx_px(sx - px);
   interval_number sy_py(sy - py);
   interval_number sz_pz(sz - pz);
   interval_number tmp_a(qx_px * ry_py);
   interval_number tmp_b(qy_py * rx_px);
   interval_number m01(tmp_a - tmp_b);
   interval_number tmq_a(qx_px * rz_pz);
   interval_number tmq_b(qz_pz * rx_px);
   interval_number m02(tmq_a - tmq_b);
   interval_number tmr_a(qy_py * rz_pz);
   interval_number tmr_b(qz_pz * ry_py);
   interval_number m12(tmr_a - tmr_b);
   interval_number mt1(m01 * sz_pz);
   interval_number mt2(m02 * sy_py);
   interval_number mt3(m12 * sx_px);
   interval_number mtt(mt1 - mt2);
   interval_number m012(mtt + mt3);
   interval_number qx_px1(qx_px - dx);
   interval_number qy_py1(qy_py - dy);
   interval_number rx_px1(rx_px - dx);
   interval_number ry_py1(ry_py - dy);
   interval_number rz_pz1(rz_pz - dz);
   interval_number qz_pz1(qz_pz - dz);
   interval_number sx_px1(sx_px - dx);
   interval_number sy_py1(sy_py - dy);
   interval_number sz_pz1(sz_pz - dz);
   interval_number tmp_a1(qx_px1 * ry_py1);
   interval_number tmp_b1(qy_py1 * rx_px1);
   interval_number m011(tmp_a1 - tmp_b1);
   interval_number tmq_a1(qx_px1 * rz_pz1);
   interval_number tmq_b1(qz_pz1 * rx_px1);
   interval_number m021(tmq_a1 - tmq_b1);
   interval_number tmr_a1(qy_py1 * rz_pz1);
   interval_number tmr_b1(qz_pz1 * ry_py1);
   interval_number m121(tmr_a1 - tmr_b1);
   interval_number mt11(m011 * sz_pz1);
   interval_number mt21(m021 * sy_py1);
   interval_number mt31(m121 * sx_px1);
   interval_number mtt1(mt11 - mt21);
   interval_number m0121(mtt1 + mt31);
   interval_number res(m012 - m0121);
   setFPUModeToRoundNEAR();

   if (!res.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return res.sign();
}

int deltaVolume_exact(double px, double py, double pz, double dx, double dy, double dz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   expansionObject o;
   double qx_px[2];
   o.two_Diff(qx, px, qx_px);
   double qy_py[2];
   o.two_Diff(qy, py, qy_py);
   double rx_px[2];
   o.two_Diff(rx, px, rx_px);
   double ry_py[2];
   o.two_Diff(ry, py, ry_py);
   double rz_pz[2];
   o.two_Diff(rz, pz, rz_pz);
   double qz_pz[2];
   o.two_Diff(qz, pz, qz_pz);
   double sx_px[2];
   o.two_Diff(sx, px, sx_px);
   double sy_py[2];
   o.two_Diff(sy, py, sy_py);
   double sz_pz[2];
   o.two_Diff(sz, pz, sz_pz);
   double tmp_a[8];
   int tmp_a_len = o.Gen_Product(2, qx_px, 2, ry_py, tmp_a);
   double tmp_b[8];
   int tmp_b_len = o.Gen_Product(2, qy_py, 2, rx_px, tmp_b);
   double m01[16];
   int m01_len = o.Gen_Diff(tmp_a_len, tmp_a, tmp_b_len, tmp_b, m01);
   double tmq_a[8];
   int tmq_a_len = o.Gen_Product(2, qx_px, 2, rz_pz, tmq_a);
   double tmq_b[8];
   int tmq_b_len = o.Gen_Product(2, qz_pz, 2, rx_px, tmq_b);
   double m02[16];
   int m02_len = o.Gen_Diff(tmq_a_len, tmq_a, tmq_b_len, tmq_b, m02);
   double tmr_a[8];
   int tmr_a_len = o.Gen_Product(2, qy_py, 2, rz_pz, tmr_a);
   double tmr_b[8];
   int tmr_b_len = o.Gen_Product(2, qz_pz, 2, ry_py, tmr_b);
   double m12[16];
   int m12_len = o.Gen_Diff(tmr_a_len, tmr_a, tmr_b_len, tmr_b, m12);
   double mt1[64];
   int mt1_len = o.Gen_Product(m01_len, m01, 2, sz_pz, mt1);
   double mt2[64];
   int mt2_len = o.Gen_Product(m02_len, m02, 2, sy_py, mt2);
   double mt3[64];
   int mt3_len = o.Gen_Product(m12_len, m12, 2, sx_px, mt3);
   double mtt[128];
   int mtt_len = o.Gen_Diff(mt1_len, mt1, mt2_len, mt2, mtt);
   double m012_p[128], *m012 = m012_p;
   int m012_len = o.Gen_Sum_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 128);
   double qx_px1[3];
   double qy_py1[3];
   double rx_px1[3];
   double ry_py1[3];
   double rz_pz1[3];
   double qz_pz1[3];
   double sx_px1[3];
   double sy_py1[3];
   double sz_pz1[3];
   double tmp_a1[18];
   int tmp_a1_len = o.Gen_Product(3, qx_px1, 3, ry_py1, tmp_a1);
   double tmp_b1[18];
   int tmp_b1_len = o.Gen_Product(3, qy_py1, 3, rx_px1, tmp_b1);
   double m011[36];
   int m011_len = o.Gen_Diff(tmp_a1_len, tmp_a1, tmp_b1_len, tmp_b1, m011);
   double tmq_a1[18];
   int tmq_a1_len = o.Gen_Product(3, qx_px1, 3, rz_pz1, tmq_a1);
   double tmq_b1[18];
   int tmq_b1_len = o.Gen_Product(3, qz_pz1, 3, rx_px1, tmq_b1);
   double m021[36];
   int m021_len = o.Gen_Diff(tmq_a1_len, tmq_a1, tmq_b1_len, tmq_b1, m021);
   double tmr_a1[18];
   int tmr_a1_len = o.Gen_Product(3, qy_py1, 3, rz_pz1, tmr_a1);
   double tmr_b1[18];
   int tmr_b1_len = o.Gen_Product(3, qz_pz1, 3, ry_py1, tmr_b1);
   double m121[36];
   int m121_len = o.Gen_Diff(tmr_a1_len, tmr_a1, tmr_b1_len, tmr_b1, m121);
   double mt11_p[128], *mt11 = mt11_p;
   int mt11_len = o.Gen_Product_With_PreAlloc(m011_len, m011, 3, sz_pz1, &mt11, 128);
   double mt21_p[128], *mt21 = mt21_p;
   int mt21_len = o.Gen_Product_With_PreAlloc(m021_len, m021, 3, sy_py1, &mt21, 128);
   double mt31_p[128], *mt31 = mt31_p;
   int mt31_len = o.Gen_Product_With_PreAlloc(m121_len, m121, 3, sx_px1, &mt31, 128);
   double mtt1_p[128], *mtt1 = mtt1_p;
   int mtt1_len = o.Gen_Diff_With_PreAlloc(mt11_len, mt11, mt21_len, mt21, &mtt1, 128);
   double m0121_p[128], *m0121 = m0121_p;
   int m0121_len = o.Gen_Sum_With_PreAlloc(mtt1_len, mtt1, mt31_len, mt31, &m0121, 128);
   double res_p[128], *res = res_p;
   int res_len = o.Gen_Diff_With_PreAlloc(m012_len, m012, m0121_len, m0121, &res, 128);

   double return_value = res[res_len - 1];
   if (res_p != res) free(res);
   if (m0121_p != m0121) free(m0121);
   if (mtt1_p != mtt1) free(mtt1);
   if (mt31_p != mt31) free(mt31);
   if (mt21_p != mt21) free(mt21);
   if (mt11_p != mt11) free(mt11);
   if (m012_p != m012) free(m012);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

int deltaVolume(double px, double py, double pz, double dx, double dy, double dz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   int ret;
   ret = deltaVolume_filtered(px, py, pz, dx, dy, dz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = deltaVolume_interval(px, py, pz, dx, dy, dz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return deltaVolume_exact(px, py, pz, dx, dy, dz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
}

bool rayIntersectsPlane(const explicitPoint3D& src, const explicitPoint3D& dir,
	const explicitPoint3D& v1, const explicitPoint3D& v2, const explicitPoint3D& v3)
{
    int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);

    int o1 = orient3d(src.X(), src.Y(), src.Z(), v1.X(), v1.Y(), v1.Z(), v2.X(), v2.Y(), v2.Z(), v3.X(), v3.Y(), v3.Z());
    return (o1 == 0 || o1 == deltaVolume(src.X(), src.Y(), src.Z(), dir.X(), dir.Y(), dir.Z(), v1.X(), v1.Y(), v1.Z(), v2.X(), v2.Y(), v2.Z(), v3.X(), v3.Y(), v3.Z()));
}

void simpleTest()
{
    explicitPoint3D s(0, 0, 0), d1(1, 0, 0), d2(-1, 0, 0);
    explicitPoint3D v1(2, -1, -1), v2(2, 1, -1), v3(2, 1, 1);

    std::cout << rayIntersectsPlane(s, d1, v1, v2, v3) << std::endl; // Should output 1
    std::cout << rayIntersectsPlane(s, d2, v1, v2, v3) << std::endl; // Should output 0
}
