from enum import Enum

class OutputFormat(str, Enum):
    csv = "csv"
    table = "table"
    excel = "excel"

