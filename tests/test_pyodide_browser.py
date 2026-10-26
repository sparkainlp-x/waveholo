import http.server
import os
import socketserver
import threading
import unittest

try:
    from playwright.sync_api import sync_playwright
    PLAYWRIGHT_AVAILABLE = True
except ImportError:
    PLAYWRIGHT_AVAILABLE = False


@unittest.skipUnless(PLAYWRIGHT_AVAILABLE, "Playwright not installed")
class TestPyodideBrowserSandbox(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        docs_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "docs", "sandbox")
        handler = lambda *a, **kw: http.server.SimpleHTTPRequestHandler(*a, directory=docs_dir, **kw)
        cls.httpd = socketserver.TCPServer(("127.0.0.1", 0), handler)
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def test_headless_pyodide_simulation(self):
        with sync_playwright() as p:
            try:
                browser = p.chromium.launch(headless=True)
            except Exception as exc:
                self.skipTest(f"Playwright browser binary unavailable: {exc}")
            page = browser.new_page()
            page.goto(f"http://127.0.0.1:{self.port}/index.html")
            page.wait_for_selector("#run-btn:not([disabled])", timeout=45000)
            page.click("#run-btn")
            page.wait_for_function("document.getElementById('console-out').innerText.includes('WaveHolo Measured Focus Plane')", timeout=30000)
            text = page.locator("#console-out").inner_text()
            self.assertIn("WaveHolo Measured Focus Plane", text)
            self.assertIn("-5.000 mm", text)
            browser.close()


if __name__ == "__main__":
    unittest.main()
