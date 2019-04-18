"""
CardioTRaNs 2019
Author: Rohit Suratekar

Kivy GUI functions
"""

from SecretColors import palette
from kivy.app import App
from kivy.lang import Builder
from kivy.properties import BooleanProperty
from kivy.uix.behaviors import FocusBehavior
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.recycleboxlayout import RecycleBoxLayout
from kivy.uix.recycleview import RecycleView
from kivy.uix.recycleview.layout import LayoutSelectionBehavior
from kivy.uix.recycleview.views import RecycleDataViewBehavior

from local.operations import LocalDatabase

p = palette.Palette(color_mode="rgba")

Builder.load_file('visualization/design.ky')


class SelectableRecycleBoxLayout(FocusBehavior,
                                 LayoutSelectionBehavior,
                                 RecycleBoxLayout):
    """
    Adds selection and focus behaviour to the view.
    """


class SelectableLabel(RecycleDataViewBehavior, Label):
    """
    Add selection support to the Label
    """
    index = None
    selected = BooleanProperty(False)
    selectable = BooleanProperty(True)
    mock_select = False
    row_color_normal = p.gray(shade=40)
    row_color_selected = p.black()

    def refresh_view_attrs(self, rv, index, data):
        """
        Catch and handle the view changes
        """
        self.index = index
        return super(SelectableLabel, self).refresh_view_attrs(
            rv, index, data)

    def on_touch_down(self, touch):
        """
        Add selection on touch down
        """
        if super(SelectableLabel, self).on_touch_down(touch):
            return True
        if self.collide_point(*touch.pos) and self.selectable:
            return self.parent.select_with_touch(self.index, touch)

    def apply_selection(self, rv, index, is_selected):
        """
        Respond to the selection of items in the view.
        """
        self.selected = is_selected
        if is_selected:
            rv.item_selected(index)
            print(index)


class RV(RecycleView):
    selected_item = -1
    original_data = []
    base_color = p.red()

    def __init__(self, **kwargs):
        super(RV, self).__init__(**kwargs)
        db = LocalDatabase()
        for i in db.get_all_anatomy_items():
            self.data.append({
                "text": "{} : {}".format(i.id, i.name),
                "is_selected": False,
                "b_color": self.base_color}
            )
        self.original_data = [x for x in self.data]

    def item_selected(self, item: int):
        if self.selected_item != item:
            self.selected_item = item


class Test(BoxLayout):
    box_background = p.gray_cool(shade=10)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def on_search(text: str, rv: RV):
        rv.selected_item = -1
        new_data = []
        if len(text.strip()) > 0:
            for r in rv.original_data:
                if text.strip().lower() in r["text"]:
                    new_data.append(r)
        else:
            new_data.extend(rv.original_data)

        rv.data = new_data

    @staticmethod
    def clear_search(text_input):
        text_input.text = ""


class TestApp(App):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def build(self):
        return Test()


def test():
    TestApp().run()


def run():
    test()
