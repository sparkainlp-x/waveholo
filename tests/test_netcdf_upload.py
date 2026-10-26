import struct
import unittest
import numpy as np


def create_netcdf3_bytes(z_arr, r_arr, title="WaveHolo Focus Curve Simulation"):
    header = bytearray(b"CDF\x01") + struct.pack(">I", 0)
    header += struct.pack(">II", 0x0A, 1) + struct.pack(">I", 7) + b"z_sweep\x00" + struct.pack(">I", len(z_arr))
    value = title.encode("ascii")
    header += struct.pack(">II", 0x0C, 1) + struct.pack(">I", 5) + b"title\x00\x00\x00"
    header += struct.pack(">II", 2, len(value)) + value + b"\0" * ((4 - len(value) % 4) % 4)
    header += struct.pack(">II", 0x0B, 2)
    for name in (b"distance_z", b"realness"):
        pad = b"\0" * ((4 - len(name) % 4) % 4)
        header += struct.pack(">I", len(name)) + name + pad
        header += struct.pack(">IIIII", 1, 0, 0, 6, len(z_arr) * 8) + struct.pack(">I", 0)
    header += b"\0" * ((4 - len(header) % 4) % 4)
    return bytes(header + b"".join(struct.pack(">d", float(v)) for v in (*z_arr, *r_arr)))


def parse_netcdf_upload(raw_bytes, n_points=101):
    if raw_bytes[:4] != b"CDF\x01":
        raise ValueError("Invalid NetCDF format")
    expected = n_points * 16
    if len(raw_bytes) < 4 + expected:
        raise ValueError("Truncated file")
    data = raw_bytes[-expected:]
    z = np.array(struct.unpack(f">{n_points}d", data[:n_points * 8]))
    r = np.array(struct.unpack(f">{n_points}d", data[n_points * 8:]))
    i = int(np.argmin(r))
    return {"z_sweep": z, "realness": r, "best_idx": i,
            "best_z": float(z[i]), "min_realness": float(r[i])}


class TestNetCDFUploadParsing(unittest.TestCase):
    def setUp(self):
        self.sweep = np.linspace(-10.0, 0.0, 101) * 1e-3
        self.realness = 0.5 * (1 - 0.9 / (1 + ((self.sweep + 4.998e-3) / 0.2e-3) ** 2))

    def test_valid_file_decoding(self):
        result = parse_netcdf_upload(create_netcdf3_bytes(self.sweep, self.realness), len(self.sweep))
        self.assertEqual(len(result["z_sweep"]), 101)
        np.testing.assert_allclose(result["z_sweep"], self.sweep, atol=1e-15)
        np.testing.assert_allclose(result["realness"], self.realness, atol=1e-15)

    def test_focus_plane_accuracy(self):
        result = parse_netcdf_upload(create_netcdf3_bytes(self.sweep, self.realness), 101)
        self.assertAlmostEqual(result["best_z"] * 1e3, -5.0, delta=0.1)
        self.assertLess(result["min_realness"], 0.1)

    def test_bad_magic_rejected(self):
        with self.assertRaises(ValueError):
            parse_netcdf_upload(b"not-a-netcdf")

    def test_truncated_file_rejected(self):
        raw = create_netcdf3_bytes(self.sweep, self.realness)
        with self.assertRaises(ValueError):
            parse_netcdf_upload(raw[:100], 101)


if __name__ == "__main__":
    unittest.main()
