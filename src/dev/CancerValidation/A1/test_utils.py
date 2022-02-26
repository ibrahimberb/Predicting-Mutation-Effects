from unittest import TestCase
from utils import (
    CancerValidation,
    EntryNotFoundError,
    ClassLabels,
    unzip_res_range
)


class TestCancerValidation(TestCase):
    H_SAPIENS_INTERFACES_HQ_PATH = "../data/H_sapiens_interfacesHQ.txt"

    def setUp(self) -> None:
        self.cancer_validation = CancerValidation(
            interfaces_data_path=self.H_SAPIENS_INTERFACES_HQ_PATH
        )

    # def test_check(self):
    #     self.fail()

    # def test__get_entry(self):
    #     self.fail()

    def test_res(self):
        _, res = self.cancer_validation._get_entry(
            protein="Q71DI3", interactor="P62805"
        )

        self.assertEqual(
            res,
            unzip_res_range(
                "[59-64,66-67,70-71,74-75,79-80,82-86,88-89,92-93,96-98,100-103,105-107,109,118-122,"
                "125-126,129,131-132,134]"
            )
        )

    def test_entry_not_found(self):
        with self.assertRaises(EntryNotFoundError):
            self.cancer_validation.check(
                protein="Q96GD4", mutation="I128T", interactor="Q9ULW0"
            )

    def test_residue_not_found(self):
        self.assertEqual(
            self.cancer_validation.check(
                protein="P03951", mutation="R522K", interactor="P14210"
            ),
            ClassLabels.NON_DISRUPTIVE,
        )
